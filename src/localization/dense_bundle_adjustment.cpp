/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

/*
Cost functions for localization based on photometric data and mutual information
*/

#include "localization/local_cost_functions.h"

#include "std.h"
#include "eigen.h"
#include "ocv.h" //TODO replace here Mat32f by a pointer
#include "ceres.h"
#include "io.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "projection/jacobian.h"
#include "reconstruction/triangulator.h"

DenseBaCostFunction::DenseBaCostFunction(const ICamera * camera, const Transf & xiBaseCam,
            const PhotometricPack & dataPack,
            const Mat32f & img1, const Mat32f & img2, double scale) :
            _camera(camera->clone()),
            _dataPack(dataPack),
            _xiBaseCam(xiBaseCam),
            _imageGrid1(img1.cols, img1.rows, (float*)(img1.data)),
            _imageGrid2(img2.cols, img2.rows, (float*)(img2.data)),
            _invScale(1. / scale),
            MARGIN_SIZE(50. / scale)
    {
        mutable_parameter_block_sizes()->clear();
        mutable_parameter_block_sizes()->push_back(6);
        set_num_residuals(_dataPack.cloud.size());
    }
    
void DenseBaCostFunction::lossFunction(const double x, double & rho, double & drhodx) const
{
    double s = sign(x);
    const double e = exp(-abs(x) / LOSS_FACTOR);
    rho = s * LOSS_FACTOR * (1. - e);
    drhodx = e;

    return;
}

/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera
*/

//TODO derive the class from Photometric cost 
double DenseBaCostFunction::getUMapgin(const double & u) const
{
    double uScale = u * _invScale;
    if (uScale < MARGIN_SIZE) 
    {
        return uScale - MARGIN_SIZE;
    }
    else if (uScale > _imageGrid.uMax - MARGIN_SIZE - 1)
    {
        return uScale - _imageGrid.uMax + MARGIN_SIZE + 1;
    } 
    else return 0;
}

double DenseBaCostFunction::getVMapgin(const double & v) const
{
    double vScale = v * _invScale;
    if (vScale < MARGIN_SIZE) 
    {
        return vScale - MARGIN_SIZE;
    }
    else if (vScale > _imageGrid.vMax - MARGIN_SIZE - 1)
    {
        return vScale - _imageGrid.vMax + MARGIN_SIZE + 1;
    } 
    else return 0;
}

bool DenseBaCostFunction::Evaluate(double const * const * parameters,
        double * residual, double ** jacobian) const
{
    const int POINT_NUMBER = _dataPack.cloud.size();
    
    for (int imgIdx = 0; imgIdx < 2; imgIdx++)
    {
        Transf xiBase(parameters[imgIdx]);
        Transf xiCam = xiBase.compose(_xiBaseCam);
        // point cloud in frame 2
        vector<Vector3d> transformedPoints = _dataPack.cloud;
        for (int i = 0; i < POINT_NUMBER; i++)
        {
            transformedPoints[i] *= parameters[2][i];
        }
        xiCam.inverseTransform(transformedPoints, transformedPoints);
        
        // init the image interpolation
        ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(imgIdx ? _imageGrid2 : _imageGrid1);
        
        bool computeJac = (jacobian != NULL);
        const double FADE = 0.01 * MARGIN_SIZE * MARGIN_SIZE;
        if (computeJac)
        {
            // L_uTheta
            CameraJacobian jacobianCalculator(_camera, xiBase, _xiBaseCam);
            for (int i = 0; i < POINT_NUMBER; i++)
            {
                Vector2d pt;
    //            bool projRes = 
                if (not _camera->projectPoint(transformedPoints[i], pt)) 
                {
                    residual[i] = 0;
                    fill(jacobian[0] + i*6, jacobian[0] + i*6 + 6, 0.);
                    continue;
                }
                
//                const double uMarg = getUMapgin(pt[0]);
//                const double vMarg = getVMapgin(pt[1]);
                
                double uMarg = 0;
                double vMarg = 0;
                                  
                double f;
                // image interpolation and gradient
                Covector2d grad;
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale,
                        &f, &grad[1], &grad[0]);
                grad *= _invScale;  // normalize according to the scale
                
                residual[i] = (f - _dataPack.valVec[i]);
                double drhoderr;
                lossFunction(residual[i], residual[i], drhoderr);
                
                
                
                if (uMarg == 0 and vMarg == 0)
                {
                    Covector6d dfdxi;
                    jacobianCalculator.dfdxi(transformedPoints[i], grad, dfdxi.data());
                    dfdxi *= drhoderr;
                    copy(dfdxi.data(), dfdxi.data() + 6, jacobian[0] + i*6);
                }
                else
                {
                    Covector6d dudxi, dvdxi;
                    jacobianCalculator.dpdxi(transformedPoints[i], dudxi.data(), dvdxi.data());
                    Covector6d drhodxi = drhoderr * (grad[0] * dudxi + grad[1] * dvdxi); 
                
                    //fade-away margins
                    const double phi = FADE / (FADE + uMarg * uMarg + vMarg * vMarg);
                    const double K = -2 * phi * phi / FADE * _invScale * _invScale;
                    const double dphidu = K * uMarg;
                    const double dphidv = K * vMarg;
                    
                    Map<Covector6d>(jacobian[0] + i*6) = (drhodxi * phi + 
                                                        residual[i] * (dphidu * dudxi + dphidv * dvdxi))*0;
                    residual[i] *= phi * 0;
                }
            }
        }
        else
        {
            for (int i = 0; i < POINT_NUMBER; i++)
            {
                Vector2d pt;
                if (not _camera->projectPoint(transformedPoints[i], pt)) 
                {
                    residual[i] = 0.;
                    continue;
                }
                
                double f;
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale, &f);
                residual[i] = (f - _dataPack.valVec[i]);
                double k;
                lossFunction(residual[i], residual[i], k);
                const double uMarg = getUMapgin(pt[0]);
                const double vMarg = getVMapgin(pt[1]);
                if (uMarg != 0 or vMarg != 0)
                {
                    
                    const double phi = FADE / (FADE + uMarg * uMarg + vMarg * vMarg);
                    residual[i] *= phi * 0;
                }
                
            }
        }
    }
    return true;
}


