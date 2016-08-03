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
#include "ocv.h"

#include "camera/generic_camera.h"
#include "utils/scale_parameters.h"

const double DEFAULT_DEPTH = 1;
const double MIN_DEPTH = 0.1;
const double DEFAULT_SIGMA_DEPTH = 100;

class DepthMap : private ScaleParameters
{
public:
    DepthMap() : cameraPtr(NULL) {}
    
    //copy constructor
    DepthMap(const DepthMap & depth) :
            ScaleParameters(depth),
            cameraPtr(depth.cameraPtr->clone()),
            valVec(depth.valVec),
            sigmaVec(depth.sigmaVec) {}
    
    //basic constructor        
    DepthMap(const ICamera * camera, const ScaleParameters & params) :
            ScaleParameters(params),
            cameraPtr(camera->clone()),
            valVec(xMax*yMax, DEFAULT_DEPTH),
            sigmaVec(xMax*yMax, DEFAULT_SIGMA_DEPTH)  {}

    virtual ~DepthMap() 
    {
        delete cameraPtr;
        cameraPtr = NULL;
    }
    
    DepthMap & operator = (const DepthMap & other)
    {
        if (this != &other)
        {
            delete cameraPtr;
            cameraPtr = other.cameraPtr->clone();
            ScaleParameters::operator = (other);
            valVec = other.valVec;
            sigmaVec = other.sigmaVec;
        }
        return *this;
    }
    
    void setDefault()
    {
        setTo(DEFAULT_DEPTH, DEFAULT_SIGMA_DEPTH);
    }
    
    void setTo(double val, double sigmaVal)
    {
        fill(valVec.begin(), valVec.end(), val);
        fill(sigmaVec.begin(), sigmaVec.end(), sigmaVal);
    }
    
    void applyMask(const Mat8u & mask);
    
    //check the limits
    bool isValid(int x, int y) const;
    
    // nearest neighbor interpolation
    const double & nearest(int u, int v) const;
    const double & nearest(Vector2d pt) const;
    double & nearest(Vector2d pt);
    
    // nearest neighbor interpolation for the uncertainty
    const double & nearestSigma(int u, int v) const;
    const double & nearestSigma(Vector2d pt) const;
    double & nearestSigma(Vector2d pt);
    
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
    
    Vector2dVec getPointVec(const std::vector<int> idxVec) const;
    Vector2dVec getPointVec() const;
    
    void reconstructUncertainty(std::vector<int> & idxVec, 
            Vector3dVec & minDistVec,
            Vector3dVec & maxDistVec) const;
            
    void reconstruct(std::vector<int> & idxVec, Vector3dVec & result) const;
    
    // idxVec corresponds to indices of points in queryPointVec
    void reconstruct(const Vector2dVec & queryPointVec,
            std::vector<int> & idxVec, Vector3dVec & result) const;
    
    //TODO make it bool and make it return a mask
    void project(const Vector3dVec & pointVec, Vector2dVec & result) const;
    
    //TODO make a wrap method
    void toMat(Mat32f & out) const;
    
    int getWidth() const { return xMax; }
    int getHeight() const { return yMax; }
    
    static DepthMap generatePlane(const ICamera * camera, const ScaleParameters & params, 
            Transformation<double> TcameraPlane, const Vector3dVec & polygonVec);
    
private:
    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    double foo; //to return in case of out-of-range 
    
    ICamera * cameraPtr;
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

