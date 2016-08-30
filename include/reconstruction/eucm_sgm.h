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
Semi-global block matching algorithm for non-rectified images
NOTE:
(u, v) is an image point 
(x, y) is a depth map point 
*/

#pragma once

#include "std.h"
#include "ocv.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "camera/eucm.h"

#include "utils/scale_parameters.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/epipoles.h"
#include "reconstruction/eucm_stereo.h"

struct SGMParameters : public StereoParameters
{
    // cost parameters
    int lambdaStep = 5;
    int lambdaJump = 32;
    
    bool imageBasedCost = true;
    
    //TODO get a meaningful threshold/ method to discard
    //temporarily descriptor step must be 1
    //all non-salient points are discarded
    bool salientPoints = true;
    
    bool useUVCache = true;
};

class EnhancedSGM : private EnhancedStereo
{
public:
    
    EnhancedSGM(Transformation<double> T12, const EnhancedCamera * cam1,
            const EnhancedCamera * cam2, const SGMParameters & parameters) :
            EnhancedStereo(cam1, cam2, parameters),
            // initialize members
            params(parameters),
            epipolar(T12, cam1, cam2, 2500),
            epipoles(cam1, cam2, T12)
    { 
        setTransformation(T12);
        assert(params.dispMax % 2 == 0);
        createBuffer();
        computeReconstructed();
        computeRotated();
        computePinf();
        if (params.useUVCache) computeUVCache();
    }
    
    ~EnhancedSGM()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
    }
    
    // precompute coordinates for different disparities to speedup the computation
    void computeUVCache();
    
    // An interface function
    void computeStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depthMap);
    
    // An interface function
    void computeStereo(const Mat8u & img1, const Mat8u & img2, Mat32f & depthMat);
    
    //// EPIPOLAR GEOMETRY
    
    // computes reconstVec -- reconstruction of every pixel of the first image
    void computeReconstructed();
    
    // computes reconstRotVec -- reconstVec rotated into the second frame
    void computeRotated();
       
    // computes pinfVec -- projections of all the reconstructed points from the first image
    // onto the second image as if they were at infinity
    void computePinf();
    
    // calculate the coefficients of the polynomials for all the 
    void computeEpipolarIndices();
    
    //// DYNAMIC PROGRAMMING
    void createBuffer();
       
    // fill up the error buffer using 2*S-1 pixs along epipolar lines as local desctiprtors
    void computeCurveCost(const Mat8u & img1, const Mat8u & img2);
    
    void computeDynamicProgramming();
    
    void computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost);
    
    void reconstructDisparity();  // using the result of the dynamic programming
    
    //// MISCELLANEOUS
    
    // index of an object in a linear array corresponding to pixel [row, col] 
    int getLinearIndex(int x, int y) { return params.xMax*y + x; }
      
    CurveRasterizer<int, Polynomial2> getCurveRasteriser1(int idx);
    CurveRasterizer<int, Polynomial2> getCurveRasteriser2(int idx);
    
    // reconstruction
    //TODO move elsewhere (e.g. create a class stereo system with two cameras and a transformation)
    void computeDepth(Mat32f & distanceMat);
            
    double computeDepth(int x, int y, int h = 0);
    bool computeDepth(double & dist, double & sigma, int x, int y, int h = 0);
    
    void fillGaps(uint8_t * const data, const int step);
    
    int getHalfLength() { return min(4, max(params.scale - 1, 1)); }
    
    Mat32s & disparity() { return smallDisparity; }
private:
    const EnhancedEpipolar epipolar;
    const StereoEpipoles epipoles;
    
    std::vector<bool> maskVec;
    
    Vector2dVec pointVec1;  // the depth points on the image 1
    Vector3dVec reconstVec;  // reconstruction of every pixel by cam1
    Vector3dVec reconstRotVec;  // reconstVec rotated into the second frame
    Vector2dVec pinfVec;  // projection of reconstRotVec by cam2
    
    // discretized version
    Vector2iVec pointPxVec1;
    Vector2iVec pinfPxVec;
    
    // to be able to change the jump cost
    int jumpCost;
    
    const int DISPARITY_MARGIN = 20;
    Mat32s uCache, vCache;
    Mat8u errorBuffer;
    Mat8u costBuffer; //TODO maybe merge with salientBuffer
    Mat8u salientBuffer; 
    Mat32s tableauLeft, tableauRight; //FIXME check the type through the code
    Mat32s tableauTop, tableauBottom;
    Mat32s smallDisparity;
    Mat32s finalErrorMat;
    
    const SGMParameters params;
};

