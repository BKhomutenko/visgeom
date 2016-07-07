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
Depth container
*/
#pragma once

#include "std.h"
#include "eigen.h"

#include "camera/generic_camera.h"

const double DEFAULT_DEPTH = 1;
const double DEFAULT_SIGMA_DEPTH = 100;

class DepthMap
{
public:
    DepthMap() : cameraPtr(NULL) {}
    
    //copy constructor
    DepthMap(const DepthMap & depth) :
            cameraPtr(depth.cameraPtr->clone()),
            width(depth.width),
            height(depth.height),
            u0(depth.u0),
            v0(depth.v0),
            scale(depth.scale),
            valVec(depth.width*depth.height, DEFAULT_DEPTH),
            sigmaVec(depth.width*depth.height, DEFAULT_SIGMA_DEPTH)  {}
    
    //basic constructor        
    DepthMap(const ICamera * camera, int w, int h, double u0, double v0, double scale) :
            cameraPtr(camera->clone()), width(w), height(h), u0(u0), v0(v0), scale(scale),
            valVec(w*h, DEFAULT_DEPTH),  sigmaVec(w*h, DEFAULT_SIGMA_DEPTH)  {}

    virtual ~DepthMap() 
    {
        if (cameraPtr != NULL)
        {
            delete cameraPtr;
        }
    }
    
    // nearest neighbor interpolation
    double nearest(double u, double v);
    double nearest(Vector2d pt);
    
    // to access the elements directly
    double & at(int x, int y);
    const double & at(int x, int y) const;
    
    // to access the uncertainty directly
    double & sigma(int x, int y);
    const double & sigma(int x, int y) const;
    
    // image coordinates of depth points
    double u(int x);
    double v(int y);
    
    // image coordinates of the block corner
    int uc(int x);
    int vc(int y);
    
    // depth coordinates of image points
    int x(double u);
    int y(double v);
    
    void reconstruct(Vector3dVec & result);
    void reconstruct(const vector<int> & indexVec, Vector3dVec & result);
    void reconstruct(const Vector2dVec & pointVec, Vector3dVec & result);
    
    int getWidth() { return width; }
    int getHeight() { return height; }
    
    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    
private:
    ICamera * cameraPtr;
    int width;
    int height;
    
    double u0, v0; // image coordinates of the [0, 0] point
    double scale; // normally > 1, x = (u - u0) / ration 
};
