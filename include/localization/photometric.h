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

#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/eucm_stereo.h"

using ceres::BiCubicInterpolator;

struct Grid2D
{
    enum { DATA_DIMENSION = 1 };
    
    Grid2D(int uMax, int vMax, const float * const data) :
            data(data), uMax(uMax), vMax(vMax) {}
            
    void GetValue(int v, int u, double* val) const
    {
        if (u < 0 or u >= uMax or v < 0 or v >= vMax) 
        {
            *val = 0.;
            return;
        }
        *val = double(data[v*uMax + u]);
    }
    
    int uMax, vMax;
    const float * const data;
};

// to store the data for the photometric optimization
struct PhotometricPack
{
    vector<float> colorVec;
    vector<Vector3d> cloud;
    vector<int> idxVec;
    int scaleIdx;
};

template<template<typename> class Projector>
struct PhotometricError 
{
    PhotometricError(vector<double> & projectionParams, const PhotometricPack & dataPack,
            const Mat32f & img2)
    : _projectionParams(projectionParams), _dataPack(dataPack), _img2(img2) {}
            
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
        BiCubicInterpolator<Grid2D> imageIterpolator(imageGrid);
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
                imageIterpolator.Evaluate(pt[1], pt[0], &res);
                residual[i] = T(_dataPack.colorVec[i]) - res;
            }
            else
            {
                residual[i] = T(0.);
            }
        }
        return true;
    }
    
    vector<double> & _projectionParams;
    const PhotometricPack & _dataPack;
    const Mat32f & _img2;
};


//TODO implement multiscale
class PhotometricLocalization
{
public:
    PhotometricLocalization(const double * params1, const double * params2,
            const StereoParameters & stereoParams) : 
            cam1(stereoParams.imageWidth, stereoParams.imageHeight, params1),
            cam2(stereoParams.imageWidth, stereoParams.imageHeight, params2),
            params(stereoParams) 
            {
                params.init();
            }
            
    virtual ~PhotometricLocalization() {}
    
    bool computePose(const Mat32f & img2, Transformation<double> & T12);
    
    bool initCloud(const Mat32f & img1, const Mat32f & dist);
    
    bool wrapImage(const Mat32f & img2, Mat32f & imgOut, Transformation<double> & T12);
    
private:
    EnhancedCamera cam1, cam2;
    PhotometricPack dataPack; // to wrap the image
    StereoParameters params;
};




