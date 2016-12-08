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
#include "utils/mh_pack.h"
#include "reconstruction/stereo_misc.h"


//TODO move to mh_pack
enum ReconstructionFlags : uint32_t
{
    QUERY_POINTS = 1,
    IMAGE_VALUES = 2,
    MINMAX = 4,
    ALL_HYPOTHESES = 8,
    DEFAULT_VALUES = 16,
    SIGMA_VALUE = 32,
    QUERY_INDICES = 64,  //  stronger than QUERY_POINTS
    INDEX_MAPPING = 128    
};

class DepthMap : public ScaleParameters
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
            valVec(xMax*yMax*hMax, OUT_OF_RANGE),
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
            ScaleParameters::operator = (other);

            delete cameraPtr;
            cameraPtr = other.cameraPtr->clone();
            
            valVec = other.valVec;
            sigmaVec = other.sigmaVec;
            costVec = other.costVec;
            hMax = other.hMax;
            hStep = other.hStep;
        }
        return *this;
    }
    
    void setDefault()
    {
        setTo(DEFAULT_DEPTH, DEFAULT_SIGMA_DEPTH, DEFAULT_COST_DEPTH);
    }
    
    void setTo(const double val, 
        const double sigmaVal, 
        const double costVal = DEFAULT_COST_DEPTH)
    {
        fill(valVec.begin(), valVec.end(), val);
        fill(sigmaVec.begin(), sigmaVec.end(), sigmaVal);
        fill(costVec.begin(), costVec.end(), costVal);
    }
    
    // X is a 3D point in the camera frame
    bool pushHypothesis(const Vector3d & X, const double sigma);
    
    // (x, y) are the depth coordinates
    bool pushHypothesis(const int x, const int y, const double depth, const double sigma);

    // (u, v) are the image coordinates
    bool pushImageHypothesis(const int u, const int v, const double depth, const double sigma);
    
    bool filterPushHypothesis(const Vector3d & X, const double sigma);
    bool filterPushHypothesis(const int x, const int y, const double depth, const double sigma);
    
    void applyMask(const Mat8u & mask);
    
    //check the limits
    bool isValid(const int x, const int y, const int h = 0) const;
    bool isValid(const Vector2d pt, const int h=0) const;
    
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
    
    Vector2dVec getPointVec(const std::vector<int> & idxVec) const;
    Vector2dVec getPointVec() const;

    //TODO overload instead of default args
    vector<int> getIdxVec(const Vector2dVec & queryPointVec = Vector2dVec()) const;
    
    //TODO - Depecrated
    void reconstructUncertainty(std::vector<int> & idxVec, 
            Vector3dVec & minDistVec,
            Vector3dVec & maxDistVec) const;
            
    //TODO - Depecrated
    void reconstruct(std::vector<int> & idxVec, Vector3dVec & result) const;
    
    //TODO - Depecrated
    // idxVec corresponds to indices of points in queryPointVec
    void reconstruct(const Vector2dVec & queryPointVec,
            std::vector<int> & idxVec, Vector3dVec & result) const;

    // Reconstructs 2d points into corresponding 3d pointcloud
    // result.idxVec corresponds to indices of points in original queryPointVec
    // To use QUERY_POINT, insert the query point vector into result.idxVec. Do not
    //   use ALL_HYPOTHESES with QUERY_POINT.
    // To use IMAGE_VALUES, insert the value vector into result.valVec. This should
    //   only be used along with QUERY_POINT
    void reconstruct(MHPack & result, const uint32_t reconstFlags = 0 ) const;
    
    //TODO make it bool and make it return a mask
    void project(const Vector3dVec & pointVec, Vector2dVec & result) const;
    
    //TODO make a wrap method
    void toMat(Mat32f & out) const;
    void toInverseMat(Mat32f & out, const int layer = 0) const;
    void sigmaToMat(Mat32f & out) const;
    
    // access methods
    int getWidth() const { return xMax; }
    int getHeight() const { return yMax; }
    int getHypMax() const { return  hMax; }
    
    static DepthMap generatePlane(const ICamera * camera, const ScaleParameters & params, 
            Transformation<double> TcameraPlane, const Vector3dVec & polygonVec);
    
    void merge(const DepthMap & depth2);

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
            const Transformation<double> T12, DepthMap& output) const;

    DepthMap wrapDepth(const Transformation<double> T12,
            const ScaleParameters & scaleParams) const;

    // Filters all hypotheses to remove noise, using either a median filter or 
    // average filter, depending on the number of matching neighbour hypotheses
    void filterNoise();

    // make sure that hypotheses are like h1 h2 0 rhather than h1 0 h2
    void regularize();

    // Returns true if the two depths and sigmas are within an acceptable tolerance of each other
    static bool match(const double v1, const double s1, const double v2, const double s2);
    // Performs a filtered merge on the input depths and sigmas
    static void filter(double & v1, double & s1, const double v2, const double s2);
  
    bool empty() { return valVec.size() == 0; }
private:

    void pixelMedianFilter(const int x, const int y, const int h, DepthMap & dst);
    void pixelAverageFilter(const Vector3iVec & matches, DepthMap & dst);

    //Swap two hypotheses at the same pixel
    void swapHypotheses(const int x, const int y, const int h1, const int h2);

    //Sort all the hypotheses at a given pixel in the order of ascending cost
    void sortHypStack(const int x, const int y);

    // Rejects hypotheses above the threshold
    void costRejection(const double rejectionThreshold = DEFAULT_COST_DEPTH + 16);

    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    std::vector<double> costVec; // hypothesis cost
    int hMax; // Number of hypotheses
    int hStep; // Step to get to the next hypothesis
    
    ICamera * cameraPtr;
};
