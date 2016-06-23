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
Abstract camera class
*/
#pragma once

#include "eigen.h"
#include "geometry/geometry.h"

class ICamera
{
public:
    std::vector<double> params;
    int width, height;

    /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const = 0;

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const = 0;

    //TODO implement the projection and distortion Jacobian
    virtual bool projectionJacobian(const Vector3d & src,
            Eigen::Matrix<double, 2, 3> & Jac) const {return 0;}

    virtual void setParameters(const double * const newParams)
    {
        copy(newParams, newParams + params.size(), params.begin());
    }
    
    ICamera(int W, int H, int numParams) : width(W), height(H), params(numParams) {}

    virtual ~ICamera() {}
    
    virtual ICamera * clone() const = 0; 
    
    bool reconstructPointCloud(const Vector2dVec & src, Vector3dVec & dst) const
    {
        dst.resize(src.size());
        bool res = true;
        for (int i = 0; i < src.size(); i++)
        {
            res &= reconstructPoint(src[i], dst[i]);
        }  
        return res;
    }
    
    bool projectPointCloud(const Vector3dVec & src, Vector2dVec & dst) const
    {
        dst.resize(src.size());
        bool res = true;
        for (int i = 0; i < src.size(); i++)
        {
            res &= projectPoint(src[i], dst[i]);
        }  
        return res;
    }
};


