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
#include "utils/MHPack.h"

const double DEFAULT_DEPTH = 1;
const double MIN_DEPTH = 0.1;
const double DEFAULT_SIGMA_DEPTH = 100;
const double DEFAULT_COST_DEPTH = 100;
const double OUT_OF_RANGE = 0.0;


enum Flags : uint32_t
{
    RECONSTRUCT_QUERY_POINTS = 1,
    RECONSTRUCTION_WITH_IMAGE_VALUES = 2,
    RECONSTRUCTION_WITH_SIGMA = 4,
    MINMAX_DISTANCE_VEC_WITH_SIGN = 8,
    ADD_ALL_HYPOTHESES = 16
};

class DepthMap : private ScaleParameters
{
public:
    DepthMap() : 
            cameraPtr(NULL),
            hMax(1) {}
    
    //copy constructor
    DepthMap(const DepthMap & depth) :
            ScaleParameters(depth),
            cameraPtr(depth.cameraPtr->clone()),
            valVec(depth.valVec),
            sigmaVec(depth.sigmaVec),
            costVec(depth.costVec),
            hMax(depth.hMax),
            hStep(depth.hStep) {}

    //basic constructor for multi-hypothesis
    DepthMap(const ICamera * camera, const ScaleParameters & params, const int hMax = 1):
            ScaleParameters(params),
            cameraPtr(camera->clone()),
            valVec(xMax*yMax*hMax, DEFAULT_DEPTH),
            sigmaVec(xMax*yMax*hMax, DEFAULT_SIGMA_DEPTH),
            costVec(xMax*yMax*hMax, DEFAULT_COST_DEPTH),
            hMax(hMax),
            hStep(xMax*yMax) {}

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
            costVec = other.costVec;
            hMax = other.hMax;
        }
        return *this;
    }
    
    void setDefault()
    {
        setTo(DEFAULT_DEPTH, DEFAULT_SIGMA_DEPTH, DEFAULT_COST_DEPTH);
    }
    
    void setTo(const double val, const double sigmaVal, const double costVal = DEFAULT_COST_DEPTH)
    {
        fill(valVec.begin(), valVec.end(), val);
        fill(sigmaVec.begin(), sigmaVec.end(), sigmaVal);
        fill(costVec.begin(), costVec.end(), costVal);
    }
    
    void applyMask(const Mat8u & mask);
    
    //check the limits
    bool isValid(const int x, const int y, const int h = 0) const;
    
    // nearest neighbor interpolation
    double nearest(const int u, const int v, const int h = 0) const;
    double nearest(const Vector2d pt, const int h = 0) const;
    
    // nearest neighbor interpolation for the uncertainty
    double nearestSigma(const int u, const int v, const int h = 0) const;
    double nearestSigma(const Vector2d pt, const int h = 0) const;

    // nearest neighbor interpolation for the hypothesis cost
    double nearestCost(const int u, const int v, const int h = 0) const;
    double nearestCost(const Vector2d pt, const int h = 0) const;
    
    // to access the elements directly
    double & at(const int x, const int y, const int h = 0);
    const double & at(const int x, const int y, const int h = 0) const;
    
    // to access the elements directly
    double & at(const int idx);
    const double & at(const int idx) const;
    
    // to access the uncertainty directly
    double & sigma(const int x, const int y, const int h = 0);
    const double & sigma(const int x, const int y, const int h = 0) const;
    
    // to access the uncertainty directly
    double & sigma(const int idx);
    const double & sigma(const int idx) const;

    // to access the hypothesis cost directly
    double & cost(const int x, const int y, const int h = 0);
    const double & cost(const int x, const int y, const int h = 0) const;
    
    // to access the hypothesis cost directly
    double & cost(const int idx);
    const double & cost(const int idx) const;
    
    Vector2dVec getPointVec(const std::vector<int> idxVec) const;
    Vector2dVec getPointVec() const;
    
    //TODO - Remove this and next 2 functions, leave only the final 2 unified reconstruct functions
    void reconstructUncertainty(std::vector<int> & idxVec, 
            Vector3dVec & minDistVec,
            Vector3dVec & maxDistVec) const;
            
    void reconstruct(std::vector<int> & idxVec, Vector3dVec & result) const;
    
    // idxVec corresponds to indices of points in queryPointVec
    void reconstruct(const Vector2dVec & queryPointVec,
            std::vector<int> & idxVec, Vector3dVec & result) const;

    void reconstructUncertainty(MHPack & result, const Flags flags = (Flags)(MINMAX_DISTANCE_VEC_WITH_SIGN | ADD_ALL_HYPOTHESES) ) const;

    void reconstruct(MHPack & result, const Flags flags = (Flags)(RECONSTRUCTION_WITH_SIGMA | ADD_ALL_HYPOTHESES) ) const;
    
    //TODO make it bool and make it return a mask
    void project(const Vector3dVec & pointVec, Vector2dVec & result) const;
    
    //TODO make a wrap method
    void toMat(Mat32f & out) const;
    
    int getWidth() const { return xMax; }
    int getHeight() const { return yMax; }
    
    static DepthMap generatePlane(const ICamera * camera, const ScaleParameters & params, 
            Transformation<double> TcameraPlane, const Vector3dVec & polygonVec);
    
    void merge(const DepthMap & depth2, const Transformation<double> T12);
private:
    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    std::vector<double> costVec; // hypothesis cost
    int hMax; // Number of hypotheses
    int hStep; // Step to get to the next hypothesis
    
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

