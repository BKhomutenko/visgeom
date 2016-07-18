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
Relative camera pose estimation based on photometric error and depth map
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"

#include <ceres/ceres.h>
#include <ceres/cubic_interpolation.h>

#include "reconstruction/depth_map.h"
#include "localization/scale_space.h"
#include "geometry/geometry.h"
#include "camera/generic_camera.h"
#include "camera/eucm.h"

struct Grid2D
{
    enum { DATA_DIMENSION = 1 };
    
    Grid2D(int uMax, int vMax, const float * const data) :
            data(data), uMax(uMax), vMax(vMax) {}
            
    void GetValue(int v, int u, double* val) const
    {
        if (u < 0 or u >= uMax or v < 0 or v >= vMax) *val = 0.;
        else *val = double(data[v*uMax + u]);
        
    }
    
    int uMax, vMax;
    const float * const data;
};


// to store the data for the photometric optimization
struct PhotometricPack
{
    vector<float> colorVec;
    vector<Vector3d> cloud;
    vector<Vector3d> gradientVec;
    vector<double> gradValVec;
    vector<int> idxVec;
    int scaleIdx;
};


/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera

Jacobian matrix of wrapped image brightness with respect
to xi -- the transformation from camera frame 2 to camera 
frame 1 transformation.
Transformation is parametrized using uTheta angle-axis

J = grad(img) * dp/dX * [ -I  hat(X) ] * L_uTheta 

- grad(img) is computed with the interpolation

- dp/dX is provided by the camera

- [ -I  hat(X) ] represents the relation between motion of a spatial point 
and camera's kinematic screw. X and [  v  ] must be expressed in frame 2
                                    [omega]
dX/dt = [ -I  hat(X) ] * [  v  ]
                         [omega]

- L_uTheta is a mapping between dxi/dt and [  v  ] = V in frame 2
                                           [omega]     
L_uTheta =  [ R21      0        ]
            [  0   R21*M_uTheta ]
            
V = L_uTheta * dxi/dt
*/
struct PhotometricCostFunction : ceres::CostFunction
{

    PhotometricCostFunction(const ICamera * camera, const PhotometricPack & dataPack,
            const Mat32f & img2, double scale) :
            _camera(camera->clone()),
            _dataPack(dataPack),
            _img2(img2),
            _scale(scale) 
    {
        mutable_parameter_block_sizes()->clear();
        mutable_parameter_block_sizes()->push_back(6);
        set_num_residuals(_dataPack.cloud.size());
    }
    
    virtual ~PhotometricCostFunction()  
    {
        delete _camera;
        _camera = NULL;
    }
    
    virtual bool Evaluate(double const * const * parameters,
            double * residual, double ** jacobian) const
    {
        Transformation<double> T12(parameters[0]);
        
        // point cloud in frame 2
        vector<Vector3d> transformedPoints;
        T12.inverseTransform(_dataPack.cloud, transformedPoints);
        
        // init the image interpolation
        Grid2D imageGrid(_img2.cols, _img2.rows, (float*)(_img2.data));
        ceres::BiCubicInterpolator<Grid2D> imageInterpolator(imageGrid);
        
        bool computeJac = (jacobian != NULL and jacobian[0] != NULL);
        
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2d pt;
            if (_camera->projectPoint(transformedPoints[i], pt)) 
            {
                double f;
                if (computeJac)
                {
                    // image interpolation and gradient
                    Matrix<double, 1, 2> grad;
                    imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale,
                            &f, &grad[1], &grad[0]);
                    grad /= _scale;  // normalize according to the scale
                    
                    // dp/dX
                    Matrix<double, 2, 3> projJac;
                    _camera->projectionJacobian(transformedPoints[i], projJac);
                    
                    // L_uTheta
                    Matrix3d R21 = T12.rotMatInv();
                    Matrix3d M21 = R21 * interOmegaRot(T12.rot());
                    
                    // computing the eqilibrium gradient
                    Vector3d gradEq3 = R21 * _dataPack.gradientVec[i];
                    Vector2d gradEq2 = (projJac * gradEq3).normalized() * _dataPack.gradValVec[i]; 
                    
                    // intermediate matrix
                    Matrix<double, 1, 3> dfdX = (grad + 0. * gradEq2.transpose()) * projJac;
                    
                    // fill up the corresponding jacobian row
                    Map<Matrix<double, 1, 3>> jacV(jacobian[0] + i*6);
                    Map<Matrix<double, 1, 3>> jacW(jacobian[0] + i*6 + 3);
                    jacV = -dfdX * R21; 
                    jacW = dfdX * hat(transformedPoints[i]) * M21;
                }
                else 
                {
                    imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale, &f);
                }
                
                residual[i] = f - _dataPack.colorVec[i];
            }
            else
            {
                residual[i] = 0.;
                if (computeJac)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        jacobian[0][i*6 + j] = 0.;
                    }
                }
            }
        }
        return true;
    }
    
    // TODO precompute equilibrium jacobian
    ICamera * _camera;
    const PhotometricPack & _dataPack;
    const Mat32f & _img2;   //FIXME potentially dangerous practice,
                            // one must be sure that the origian image exists
                            // during the whole solving process
    const double _scale;
};


/*
This struct is ment for the automatic differentiation by ceres-solver.
class Projector must be known at the compilation time.
*/
template<template<typename> class Projector>
struct PhotometricError 
{
    PhotometricError(const vector<double> & projectionParams, const PhotometricPack & dataPack,
            const Mat32f & img2, double scale) :
            _projectionParams(projectionParams),
            _dataPack(dataPack),
            _img2(img2),
            _scale(scale) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> T12(params[0]);
        vector<Vector3<T>> transformedPoints;
        transformedPoints.reserve(_dataPack.cloud.size());
        for (auto & point : _dataPack.cloud)
        {
            transformedPoints.push_back(point.template cast<T>());
        }
        T12.inverseTransform(transformedPoints, transformedPoints);
        Projector<T> projector;
        Grid2D imageGrid(_img2.cols, _img2.rows, (float*)(_img2.data));
        ceres::BiCubicInterpolator<Grid2D> imageInterpolator(imageGrid);
        vector<T> projectionParamsT;
        projectionParamsT.reserve(_projectionParams.size());
        for (auto & x : _projectionParams)
        {
            projectionParamsT.push_back(T(x));
        }
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> pt;
            if (projector(projectionParamsT.data(), transformedPoints[i].data(), pt.data())) 
            {
                T res;
                imageInterpolator.Evaluate(pt[1] / T(_scale), pt[0] / T(_scale), &res);
                residual[i] = res - T(_dataPack.colorVec[i]);
            }
            else
            {
                residual[i] = T(0.);
            }
        }
        return true;
    }
    
    const vector<double> & _projectionParams;
    const PhotometricPack & _dataPack;
    const Mat32f & _img2;
    const double _scale; 
};


//TODO add assertions
class ScalePhotometric
{
public:
    ScalePhotometric(int nScales, EnhancedCamera cam2) :
            scaleSpace1(nScales),
            scaleSpace2(nScales),
            camPtr2(cam2.clone()),
            verbosity(0) {}
            
    ScalePhotometric() : camPtr2(NULL) {}
    virtual ~ScalePhotometric()
    {
        delete camPtr2;
    }
    
    void setCamera(const EnhancedCamera & cam)
    {
        delete camPtr2;
        camPtr2 = cam.clone();
    }
    
    void setNumberScales(int numScales)
    {
        scaleSpace1.setNumberScales(numScales);
        scaleSpace2.setNumberScales(numScales);
    }
    
    const DepthMap & depth() const { return depthMap; }
    DepthMap & depth() { return depthMap; }
    
    void computeBaseScaleSpace(const Mat32f & img1);
    
    PhotometricPack initPhotometricData(int scaleIdx);
    
    void computePose(const Mat32f & img2, Transformation<double> & T12);
    
    // scaleSpace2 must be initialized
    void computePose(int scaleIdx, Transformation<double> & T12);
    void computePoseAuto(int scaleIdx, Transformation<double> & T12);
    
    void setVerbosity(int newVerbosity) { verbosity = newVerbosity; }
    
private:
    BinaryScalSpace scaleSpace1;
    BinaryScalSpace scaleSpace2;
    EnhancedCamera * camPtr2;
    DepthMap depthMap;
    
    //TODO make a parameter structure
    // minimal squared norm of gradient for a pixel to be accepted
    const double GRAD_THRESH = 400;
    const double GRAD_MAX = 255;
    int verbosity;
};



