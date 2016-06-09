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

#include <vector>

#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>

#include <ceres/ceres.h>
#include <ceres/cubic_interpolation.h>

#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/eucm_stereo.h"

using std::vector;

using Eigen::Vector2d;
using Eigen::Vector3d;

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

template<template<typename> class Projector>
struct PhotometricError 
{
    PhotometricError(const vector<double> & projectionParams,
            const vector<float> & colorVec, const vector<Vector3d> & cloud,
            const cv::Mat_<float> & img2)
    : _projectionParams(projectionParams), _colorVec(colorVec), _cloud(cloud), _img2(img2) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> T12(params[0]);
        vector<Vector3<T>> transformedPoints;
        transformedPoints.reserve(_cloud.size());
        for (auto & point : _cloud)
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
                residual[i] = T(_colorVec[i]) - res;
            }
            else
            {
                residual[i] = T(0.);
            }
        }
        return true;
    }
    
    const vector<double> _projectionParams;
    const vector<float> _colorVec;
    const vector<Vector3d> _cloud;
    const cv::Mat_<float> & _img2;
};


//TODO implement multiscale
class PhotometricLocalization
{
public:
    PhotometricLocalization(int imageWidth, int imageHeight, 
            const double * params1, const double * params2,
            const StereoParameters & stereoParams) : 
            cam1(imageWidth, imageHeight, params1),
            cam2(imageWidth, imageHeight, params2),
            blockSize(stereoParams.blockSize),
            u0(stereoParams.u0 + stereoParams.disparityMax + stereoParams.blockSize),
            v0(stereoParams.v0) {}
            
    virtual ~PhotometricLocalization() {}
    
    bool computePose(const cv::Mat_<float> & img2, Transformation<double> & T12);
    
    bool initCloud(const cv::Mat_<float> & img1, const cv::Mat_<float> & dist);
    
    bool wrapImage(const cv::Mat_<float> & img2, cv::Mat_<float> & imgOut, Transformation<double> & T12);
    
private:
    EnhancedCamera cam1, cam2;
    vector<Vector3d> cloud1;
    vector<float> colorVec;
    vector<int> indexVec; // to wrap the image
    int u0, v0, blockSize;
};




