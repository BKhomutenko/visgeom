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

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"
#include "ceres.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"



// to store the data for the photometric optimization
struct PhotometricPack
{
    vector<double> valVec;
    vector<Vector3d> cloud;
    vector<int> idxVec;
    int scaleIdx;
};

/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera

- 3D points in the dataPack must be projected into the odometry base frame
- the computed transformation will correspond the the motion of the odometry frame
*/
struct PhotometricCostFunction : ceres::CostFunction
{

    PhotometricCostFunction(const ICamera * camera, const Transf & xiBaseCam,
            const PhotometricPack & dataPack,
            const Mat32f & img2, double scale);
    
    virtual ~PhotometricCostFunction()  
    {
        delete _camera;
        _camera = NULL;
    }
    
    virtual bool Evaluate(double const * const * parameters, double * residual, double ** jacobian) const;

    void lossFunction(const double x, double & rho, double & drhodx) const;
    
    double getUMapgin(const double & u) const;
    double getVMapgin(const double & v) const;
    
    ICamera * _camera;
    const Transf _xiBaseCam;
    const PhotometricPack & _dataPack;
    const Grid2D<float> _imageGrid;
    
    
//    const double _scale;
    const double _invScale;
    const double LOSS_FACTOR = 3;     // defines how quickly the impact of data points is reduced with error
    const double MARGIN_SIZE;
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
            _imageGrid(img2.cols, img2.rows, (float*)(img2.data)),
            _invScale(1. / scale) {}
            
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
        ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(_imageGrid);
        vector<T> projectionParamsT;
        projectionParamsT.reserve(_projectionParams.size());
        for (auto & x : _projectionParams)
        {
            projectionParamsT.push_back(T(x));
        }
        for (int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> pt;
            if (projector(projectionParamsT.data(), transformedPoints[i].data(), pt.data())) 
            {
                T res;
                imageInterpolator.Evaluate(pt[1] * T(_invScale), pt[0] * T(_invScale), &res);
                residual[i] = res - T(_dataPack.valVec[i]);
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
    const Grid2D<float> _imageGrid;
//    const double _scale; 
    const double _invScale;
};

// excluding the point reconstruction reduces the number of unknowns,
// might be faster TODO check it
// but statisically not optimal
struct EssentialCost : ceres::SizedCostFunction<6, 6>
{
    EssentialCost(const Vector3d x1, const Vector3d x2);

    virtual ~EssentialCost() { }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    vector<double> _aVec;
};

//5-point algorithm
//a simpler solution, more variables TODO compare performance
//TODO change from 5 to any number
struct MonoReprojectCost : ceres::SizedCostFunction<10, 6, 5>
{
    MonoReprojectCost(const ICamera * camera,
            const Vector3dVec & xVec1, const Vector2dVec & pVec2,
            const Transf xiBaseCam):
    _pVec2(pVec2), _xVec1(xVec1), _camera(camera->clone()), _xiBaseCam(xiBaseCam) 
    {
        assert(_xVec1.size() == 5);
        assert(_pVec2.size() == 5);
    }

    virtual ~MonoReprojectCost() { delete _camera; }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    const ICamera * _camera;
    const Vector3dVec _xVec1;
    const Vector2dVec _pVec2;
    Transf _xiBaseCam;
    vector<double> _aVec;
};


struct SparseReprojectCost : ceres::CostFunction
{
    SparseReprojectCost(const ICamera * camera,
            const Vector3dVec & xVec1, const Vector3dVec & xVec2,
            const Vector2dVec & pVec2, const vector<double> & sizeVec, const Transf xiBaseCam):
        _pVec2(pVec2), _xVec1(xVec1), _xVec2(xVec2), _sizeVec(sizeVec),
        _camera(camera->clone()), _xiBaseCam(xiBaseCam) 
    {
        assert(_pVec2.size() == _xVec1.size());
        assert(_pVec2.size() == _xVec2.size());
        set_num_residuals(_pVec2.size() * 2);           // projection error
        mutable_parameter_block_sizes()->push_back(6);  // 6dof transformation
    }

    virtual ~SparseReprojectCost() { delete _camera; }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    const ICamera * _camera;
    const Vector3dVec _xVec1;
    const Vector3dVec _xVec2;
    const Vector2dVec _pVec2;
    Transf _xiBaseCam;
    vector<double> _sizeVec;
};

struct OdometryPrior : ceres::SizedCostFunction<6, 6>
{
    OdometryPrior(const double errV, const double errW, const double lambdaT, const double lambdaR,
        const Transf xiOdom);

    virtual ~OdometryPrior() { }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    Transf _xiPrior;
    Matrix6d _A;
    Matrix6d _J;
};

