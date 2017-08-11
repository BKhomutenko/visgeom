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
struct DenseBaCostFunction : ceres::CostFunction
{

    DenseBaCostFunction(const ICamera * camera, const Transf & xiBaseCam,
            const PhotometricPack & dataPack,
            const Mat32f & img1, const Mat32f & img2, double scale);
    
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
    const Grid2D<float> _imageGrid1;
    const Grid2D<float> _imageGrid2;
    
//    const double _scale;
    const double _invScale;
    const double LOSS_FACTOR = 3;     // defines how quickly the impact of data points is reduced with error
    const double MARGIN_SIZE;
};

