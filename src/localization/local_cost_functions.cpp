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

PhotometricCostFunction::PhotometricCostFunction(const ICamera * camera, const Transf & xiBaseCam,
            const PhotometricPack & dataPack,
            const Mat32f & img2, double scale) :
            _camera(camera->clone()),
            _dataPack(dataPack),
            _xiBaseCam(xiBaseCam),
            _imageGrid(img2.cols, img2.rows, (float*)(img2.data)),
            _invScale(1. / scale),
            MARGIN_SIZE(50. / scale)
    {
        mutable_parameter_block_sizes()->clear();
        mutable_parameter_block_sizes()->push_back(6);
        set_num_residuals(_dataPack.cloud.size());
    }
    
void PhotometricCostFunction::lossFunction(const double x, double & rho, double & drhodx) const
{
    double s = 0.1*sign(x);
    const double arg = -abs(x) / LOSS_FACTOR;
    const double e = arg > -5 ? exp(-abs(x) / LOSS_FACTOR) : 0;
    rho = s * LOSS_FACTOR * (1. - e);
    drhodx = 0.1*e;

    return;
}

/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera
*/

double PhotometricCostFunction::getUMapgin(const double & u) const
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

double PhotometricCostFunction::getVMapgin(const double & v) const
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

bool PhotometricCostFunction::Evaluate(double const * const * parameters,
        double * residual, double ** jacobian) const
{
    const int POINT_NUMBER = _dataPack.cloud.size();
    
    Transf xiBase(parameters[0]);
    Transf xiCam = xiBase.compose(_xiBaseCam);
    // point cloud in frame 2
    vector<Vector3d> transformedPoints;
    xiCam.inverseTransform(_dataPack.cloud, transformedPoints);
    
    // init the image interpolation
    ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(_imageGrid);
    
    bool computeJac = (jacobian != NULL and jacobian[0] != NULL);
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
            
            const double uMarg = getUMapgin(pt[0]);
            const double vMarg = getVMapgin(pt[1]);
            
//            double uMarg = 0;
//            double vMarg = 0;
                              
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
//                for (int k = 0; k < 6; k++)
//                {
//                    if (abs(jacobian[0][i*6  + k]) > 1e5)
//                    {
//                        cout << "Jac is BIG : " << jacobian[0][i*6  + k] << endl;
//                        cout << "Transformed point " << transformedPoints[i].transpose() << endl;
//                        cout << "projection " << pt.transpose() << " //// " 
//                                << _imageGrid.uMax << " " << _imageGrid.vMax << endl;
//                        cout << "gradient " << grad << endl;
//                        break;
//                    }
//                }
            }
            else
            {
                Covector6d dudxi, dvdxi;
                jacobianCalculator.dpdxi(transformedPoints[i], dudxi.data(), dvdxi.data());
                Covector6d drhodxi = (drhoderr * grad[0]) * dudxi + (drhoderr * grad[1]) * dvdxi; 
            
                //fade-away margins
                const double phi = FADE / (FADE + uMarg * uMarg + vMarg * vMarg);
                const double K = -2 * phi * phi / FADE * _invScale;
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
    return true;
}




bool MonoReprojectCost::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //compute the transformation chain
    Transf xiOdom(params[0]);
    Transf xi21 = _xiBaseCam.inverseCompose(xiOdom.inverseCompose(_xiBaseCam));
    
    //compute the points in the camera frame
    Vector3dVec xVec2 =_xVec1;
    for (int i = 0; i < 5; i++)
    {
        xVec2[i] *= params[1][i];
    }
    xi21.transform(xVec2, xVec2);
    
    //compute the reprojection error
    for (int i = 0; i < 5; i++)
    {
        Vector2d modProj;
        if (_camera->projectPoint(xVec2[i], modProj)) 
        {
            Map<Vector2d> diff(residual + 2*i);
            diff = modProj - _pVec2[i];
        }
        else
        {
            residual[2*i + 1] = residual[2*i] = DOUBLE_BIG;
        }
    }
    
    if (jacobian != NULL)
    {
        //odometry jacobian
        if (jacobian[0] != NULL) 
        {
            InterJacobian jacobianCalculator(_camera, _xiBaseCam.inverse(),
                                    xiOdom, JAC_INVERTED);
            
            for (int i = 0; i < 5; i++)
            {
                jacobianCalculator.dpdxi(xVec2[i], jacobian[0] + i*12,
                                                   jacobian[0] + i*12 + 6);
            }
        }
        
        //length jacobian
        if (jacobian[1] != NULL)
        {
            Matrix3d R21 = xi21.rotMat();
            fill(jacobian[1], jacobian[1] + 50, 0);
            for (int i = 0; i < 5; i++)
            {
                Matrix23drm dpdx;
                Vector3d n2 = R21 * _xVec1[i];
                _camera->projectionJacobian(xVec2[i], dpdx.data(), dpdx.data() + 3);
                Vector2d dpdl = dpdx * n2;
                jacobian[1][i*11] = dpdl[0];
                jacobian[1][i*11 + 5] = dpdl[1];
            }
        }
    }
    return true;
}


bool SparseReprojectCost::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //compute the transformation chain
    Transf xiOdom(params[0]);
    Transf xi12 = _xiBaseCam.inverseCompose(xiOdom.compose(_xiBaseCam));
    
    //compute the points in the camera frame
    vector<double> lambdaVec(_xVec1.size());
    vector<double> jacVec;
    Triangulator triangulator(xi12);
    if (jacobian != NULL and jacobian[0] != NULL) 
    {
        jacVec.resize(_xVec1.size() * 6);
        triangulator.computeRegular(_xVec1, _xVec2, lambdaVec.data(), NULL, jacVec.data(), NULL);
    }
    else
    {
        triangulator.computeRegular(_xVec1, _xVec2, lambdaVec.data());
    }
    
    Vector3dVec xVec2 =_xVec1;
    for (int i = 0; i < xVec2.size(); i++)
    {
        xVec2[i] *= lambdaVec[i];
    }
    xi12.inverseTransform(xVec2, xVec2);
    
    //compute the reprojection error
    for (int i = 0; i < xVec2.size(); i++)
    {
        Vector2d modProj;
        if (_camera->projectPoint(xVec2[i], modProj)) 
        {
            Map<Vector2d> diff(residual + 2*i);
            diff = (modProj - _pVec2[i]) / _sizeVec[i];
        }
        else
        {
            residual[2*i + 1] = residual[2*i] = DOUBLE_BIG;
        }
    }
    
    if (jacobian != NULL)
    {
        
        if (jacobian[0] != NULL) 
        {
            //odometry jacobian
            InterJacobian jacobianCalculator(_camera, _xiBaseCam.inverse(),
                                    xiOdom, JAC_INVERTED);
            
            for (int i = 0; i < xVec2.size(); i++)
            {
                if (residual[2*i] == DOUBLE_BIG)
                {
                    fill(jacobian[0] + i * 12, jacobian[0] + (i + 1) * 12, 0);
                }
                else
                {
                    jacobianCalculator.dpdxi(xVec2[i], jacobian[0] + i*12,
                                                   jacobian[0] + i*12 + 6);
                }
            }
            
            //length jacobian
            Matrix3d R21 = xi12.rotMatInv();
            Matrix3d RcamBase = _xiBaseCam.rotMatInv();
            
            //jacobian given by triangulate is computed wrt ( v_1_2 | om_1_2 )
            
            // om_1_2 = M * r_dot
            Matrix3d M = RcamBase * interOmegaRot(xiOdom.rot());    
            
            // v_2 = v_b + om_0_b x t_b_c
            // v_1_2 = R_c_b * t_dot + Q * r_dot
            Vector3d tBaseCam1 = RcamBase * R21.transpose() * _xiBaseCam.trans();
            Matrix3d Q = -hat(tBaseCam1) * M;
            for (int i = 0; i < xVec2.size(); i++)
            {
                if (residual[2*i] == DOUBLE_BIG) continue;
                Matrix23drm dpdx;
                Vector3d n2 = R21 * _xVec1[i];
                _camera->projectionJacobian(xVec2[i], dpdx.data(), dpdx.data() + 3);
                Vector2d dpdl = dpdx * n2;
                Map<Covector3d> dudt(jacobian[0] + i*12);
                Map<Covector3d> dudr(jacobian[0] + i*12 + 3);
                Map<Covector3d> dvdt(jacobian[0] + i*12 + 6);
                Map<Covector3d> dvdr(jacobian[0] + i*12 + 9);
                
                Map<Covector3d> dldv(jacVec.data() + i*6);
                Map<Covector3d> dldw(jacVec.data() + i*6 + 3);
                Covector3d dldt = dldv * RcamBase;
                Covector3d dldr = dldw * M + dldv * Q;
                dudt += dpdl[0] * dldt;
                dudr += dpdl[0] * dldr;
                dvdt += dpdl[1] * dldt;
                dvdr += dpdl[1] * dldr;
            }
            
            for (int i = 0; i < xVec2.size(); i++)
            {
                for (double* jptr = jacobian[0] + i*12; jptr < jacobian[0] + i*12 + 6; jptr++)
                {
                    *jptr /= _sizeVec[i];
                }
            }
        }
    }
    return true;
}

OdometryPrior::OdometryPrior(const double errV, const double errW,
        const double lambdaT, const double lambdaR,
        const Transf xiOdom) : 
    _A(Matrix6d::Zero()),
    _xiPrior(xiOdom)
{
    const double delta = xiOdom.rot()(2);
    const double l = xiOdom.trans().norm();
    
    const double delta2 = delta / 2.;
    const double l2 = l / 2.;
    
    const double s = sin(delta2);
    const double c = cos(delta2);
    
    Matrixd<3, 2> dfdu;

    dfdu <<     c,     l2 * s,  
                -s,      l2 * c, 
                0,           1;

    Matrix2d Cu;
    Cu <<   errV * errV * l * l,              0,
                      0,    errW * errW * delta * delta; 
   
    Matrix3d lambdaMat = Matrix3d::Identity();
    lambdaMat(0, 0) *= lambdaT * lambdaT;
    lambdaMat(1, 1) *= lambdaT * lambdaT;
    lambdaMat(2, 2) *= lambdaR * lambdaR;
    Matrix3d Cx = dfdu * Cu * dfdu.transpose() + lambdaMat;
    Matrix3d CxInv = Cx.inverse();
    Eigen::LLT<Matrix3d> lltOfCxInv(CxInv); // compute the Cholesky decomposition of A
    Matrix3d U = lltOfCxInv.matrixU();
//    cout << U.transpose() * U << endl;
//    cout << CxInv << endl;

//    _A.topLeftCorner<2, 2>() = U.topLeftCorner<2, 2>();
//    _A.topRightCorner<2, 1>() = U.topRightCorner<2, 1>();
//    _A(2, 2) = 1. / lambdaT;
//    _A(3, 3) = 1. / lambdaR;
//    _A(4, 4) = 1. / lambdaR;
//    _A(5, 5) = U(2, 2);
    
//    //FOR THE SIMULATION 
//    //Z -- forward, Y -- rotation
//    _A(0, 0) = U(1, 1);
//    _A(0, 2) = U(0, 1);
//    _A(2, 2) = U(0, 0);
//    _A(0, 4) = U(1, 2);
//    _A(2, 4) = U(0, 2);
//    _A(1, 1) = 1. / lambdaT;
//    _A(3, 3) = 1. / lambdaR;
//    _A(4, 4) = U(2, 2);
//    _A(5, 5) = 1. / lambdaR;
    
//        FOR THE REAL DATA 
//    Y -- forward, Z -- rotation
    _A(1, 1) = U(0, 0);
    _A(0, 0) = -U(1, 1);
    _A(0, 1) = -U(0, 1);
    _A(0, 5) = -U(1, 2);
    _A(1, 5) = U(0, 2);
    _A(2, 2) = 1. / lambdaT;
    _A(3, 3) = 1. / lambdaR;
    _A(4, 4) = 1. / lambdaR;
    _A(5, 5) = U(2, 2);

    Matrix3d M = interOmegaRot(xiOdom.rot());
    Matrix3d R = xiOdom.rotMatInv();
    
    _J.topLeftCorner<3, 3>() = _A.topLeftCorner<3, 3>() * R;
    _J.topRightCorner<3, 3>() = _A.topRightCorner<3, 3>() * R * M;
    _J.bottomLeftCorner<3, 3>() = Matrix3d::Zero();
    _J.bottomRightCorner<3, 3>() = _A.bottomRightCorner<3, 3>() * R * M;
//    cout << U << endl;
//    cout << _A << endl;
//        assert(false);
}

bool OdometryPrior::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //TODO replace the difference between two transformation by InverseCompose
    Transf xi(*params);
    Transf delta = _xiPrior.inverseCompose(xi);
    Vector6d vecDelta;
    delta.toArray(vecDelta.data());
    Map<Vector6d> res(residual);
    res = _A * vecDelta;
    
    if (jacobian != NULL)
    {
        //first transformation
        if (jacobian[0] != NULL)
        {
            Map<Matrix6drm> jac(jacobian[0]);
            jac = _J;
        }
    }
    return true;    
}
