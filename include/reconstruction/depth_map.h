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
NOTE:
(u, v) is an image point 
(x, y) is a depth map point 
*/

#pragma once

#include "io.h"
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
            valVec(depth.valVec),
            sigmaVec(depth.sigmaVec) {}
    
    //basic constructor        
    DepthMap(const ICamera * camera, int w, int h, double u0, double v0, int scale) :
            cameraPtr(camera->clone()), width(w), height(h), u0(u0), v0(v0), scale(scale),
            valVec(w*h, DEFAULT_DEPTH),  sigmaVec(w*h, DEFAULT_SIGMA_DEPTH)  {}

    virtual ~DepthMap() 
    {
        delete cameraPtr;
        cameraPtr = NULL;
    }
    
    DepthMap operator = (const DepthMap & other)
    {
        cameraPtr = other.cameraPtr->clone();
        width = other.width;
        height = other.height;
        u0 = other.u0;
        v0 = other.v0;
        scale = other.scale;
        valVec = other.valVec;
        sigmaVec = other.sigmaVec;
        return *this;
    }
    
    //check the limits
    bool isValid(int x, int y) const;
    
    // nearest neighbor interpolation
    double nearest(int u, int v) const;
    double nearest(Vector2d pt) const;
    
    // nearest neighbor interpolation for the uncertainty
    double nearestSigma(int u, int v) const;
    double nearestSigma(Vector2d pt) const;
    
    // to access the elements directly
    double & at(int x, int y);
    const double & at(int x, int y) const;
    
    // to access the elements directly
    double & at(int idx);
    const double & at(int idx) const;
    
    // to access the uncertainty directly
    double & sigma(int x, int y);
    const double & sigma(int x, int y) const;
    
    // to access the uncertainty directly
    double & sigma(int idx);
    const double & sigma(int idx) const;
    
    // image coordinates of depth points
    int u(int x) const;
    int v(int y) const;

    // depth coordinates of image points
    int x(int u) const;
    int y(int v) const;
    
    void reconstruct(Vector3dVec & result) const;
    void reconstruct(const vector<int> & indexVec, Vector3dVec & result) const;
    void reconstruct(const Vector2dVec & pointVec, Vector3dVec & result) const;
    void project(const Vector3dVec & pointVec, Vector2dVec & result) const;
    
    int getWidth() const { return width; }
    int getHeight() const { return height; }
    
private:
    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    
    ICamera * cameraPtr;
    int width;
    int height;
    
    int u0, v0; // image coordinates of the [0, 0] point
    int scale; // normally > 1, x = (u - u0) / ration 
};


class DepthReprojector
{

public:
	DepthReprojector() {}
    
    /*
    Takes the depthmap from the first image, reconstructs the pointcloud
    in the frame of the second image, then reprojects into second image
    frame. Reconstructs the cloud at the image points specified by the 
    first pointcloud, and sends this back to the first image. This 
    pointcloud is transformed into the frame of the first image, and each
    point in this cloud is projected onto the line of it's original point.
    This new depthmap is the reprojected depthmap.
    */
	void wrapDepth(const DepthMap& dMap1, const DepthMap& dMap2,
	        const Transformation<double> T12, DepthMap& output);
private:
	//null
};

