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

#include "std.h"
#include "eigen.h"
#include "geometry/geometry.h"

//TODO replace eigen vectors by double*

class ICamera
{
public:
    
    int width, height;

    /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const = 0;

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const = 0;
    
    virtual double getCenterU() { return width / 2; }
    
    virtual double getCenterV() { return height / 2; }
    
    //TODO implement the projection and distortion Jacobian
    virtual bool projectionJacobian(const Vector3d & src,
            double * dudx, double * dvdx) const {return false;}
    
    //MEMORY IS SUPPOSED TO BE ALLOCATED
    virtual bool intrinsicJacobian(const Vector3d & src,
            double * dudalpha, double * dvdalpha) const {return false;}
            
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
    
    bool reconstructPointCloud(const Vector2dVec & src,
            Vector3dVec & dst, std::vector<bool> & maskVec) const
    {
        dst.resize(src.size());
        maskVec.resize(src.size());
        bool res = true;
        for (int i = 0; i < src.size(); i++)
        {
            maskVec[i] = reconstructPoint(src[i], dst[i]);
            res &= maskVec[i];
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
    
    bool projectPointCloud(const Vector3dVec & src,
            Vector2dVec & dst, std::vector<bool> & maskVec) const
    {
        dst.resize(src.size());
        maskVec.resize(src.size());
        bool res = true;
        for (int i = 0; i < src.size(); i++)
        {
            maskVec[i] = projectPoint(src[i], dst[i]);
            res &= maskVec[i];
        }  
        return res;
    }
    
    const double * getParams() const { return params.data(); }
    
    int numParams() const { return params.size(); }
    
    virtual double lowerBound(int idx) const { return 0; }
    virtual double upperBound(int idx) const { return 1e4; }
    
protected:
    std::vector<double> params;
};


