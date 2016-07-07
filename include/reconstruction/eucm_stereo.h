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
*/

#pragma once
//STL
#include <algorithm>

#include "ocv.h"
#include "eigen.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "curve_rasterizer.h"
#include "depth_map.h"

using EpipolarRasterizer = CurveRasterizer<Polynomial2>;

struct StereoParameters
{
    // basic parameters
    int dispMax = 48; // maximum disparity
    int scale = 3;
    int uMargin = 0, vMargin = 0;  // RoI left upper corner
    int width = -1, height = -1;  // RoI size
    int lambdaStep = 5;
    int lambdaJump = 32;
    int imageWidth = 0, imageHeight = 0;
    
    // precomputed parameters
    int u0, v0;
    int uMax, vMax;
    int smallWidth, smallHeight;
    int halfBlockSize;
    
    // to be called before using
    void init()
    {
        u0 = uMargin + dispMax + scale; 
        v0 = vMargin;
        
        if (width > 0) uMax = u0 + width;
        else uMax = imageWidth - uMargin - scale;
        
        if (height > 0) vMax = v0 + height;
        else vMax = imageHeight - vMargin - scale;
        
        smallWidth = uSmall(uMax) + 1;
        smallHeight = vSmall(vMax) + 1;
        
        halfBlockSize =  scale / 2; 
    }
    
    // from image to small disparity coordiante transform
    int uSmall(int u) { return (u - u0) / scale; }
    int vSmall(int v) { return (v - v0) / scale; }
    
    // from small disparity to image coordiante transform    
    int uBig(int u) { return (u + 0.5) * scale - 0.5 + u0; }
    int vBig(int v) { return (v + 0.5) * scale - 0.5 + v0; }
    
    // from small disparity to image block corner transform 
    int uCorner(int u) { return u * scale + u0; }
    int vCorner(int v) { return v * scale + v0; }
};

class EnhancedStereo
{
public:
    EnhancedStereo(Transformation<double> T12,
            const double * params1, const double * params2, const StereoParameters & stereoParams)
            : Transform12(T12), 
            cam1(stereoParams.imageWidth, stereoParams.imageHeight, params1),
            cam2(stereoParams.imageWidth, stereoParams.imageHeight, params2),
            params(stereoParams)
    { 
        params.init();
        init(); 
    }

    void setTransformation(Transformation<double> T12) 
    { 
        Transform12 = T12;
        initAfterTransformation();
    } 
    
    void init()
    {
        createBuffer();
        computeReconstructed();
        initAfterTransformation();
    }
    
    // Only data invalidated after the transformation change are recomputed
    void initAfterTransformation()
    {
        computeEpipole();
        computeRotated();
        computePinf();
        computeEpipolarCurves();
    }
    
    //// EPIPOLAR GEOMETRY
    
    // computes reconstVec -- reconstruction of every pixel of the first image
    void computeReconstructed();
    
    // computes reconstRotVec -- reconstVec rotated into the second frame
    void computeRotated();
    
    // f2(t21) -- projection of the first camera's projection center
    void computeEpipole();
    
    // computes pinfVec -- projections of all the reconstructed points from the first image
    // onto the second image as if they were at infinity
    void computePinf();
    
    // calculate the coefficients of the polynomials for all the 
    void computeEpipolarCurves();
    
    // to visualize the epipolar lines
    void traceEpipolarLine(cv::Point pt, Mat8u & out);
    
    
    //// DYNAMIC PROGRAMMING
    
    // fills up the output with photometric errors between the val = I1(pi) and 
    // the values from I2 on the epipolar curve
    void comuteStereo(const Mat8u & img1, 
            const Mat8u & img2,
            Mat8u & disparity);
    
    void comuteStereo(const Mat8u & img1, 
            const Mat8u & img2,
            DepthMap & disparity);
            
    void createBuffer();
    
    void computeCost(const Mat8u & img1, const Mat8u & img2);
    
    void computeDynamicProgramming();
    
    void computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost);
    
    void reconstructDisparity();  // using the result of the dynamic programming
    
    // TODO implement
    void upsampleDisparity(const Mat8u & img1, Mat8u & disparity);
    
    //// MISCELLANEOUS
    
    // index of an object in a linear array corresponding to pixel [row, col] 
    //TODO change this to the disparity-based indexing
    int getLinearIdx(int u, int v) { return params.smallWidth*v + u; }
      
    // reconstruction
    Vector3d triangulate(int x1, int y1, int x2, int y2);
    void computeDistance(Mat32f & distance);
    void computeDistance(DepthMap & disparity);
    void generatePlane(Transformation<double> TcameraPlane, 
            Mat32f & distance, const Vector3dVec & polygonVec);
    
private:
    Transformation<double> Transform12;  // pose of the first to the second camera
    EnhancedCamera cam1, cam2;
   
    //TODO replace by disparity map points
    Vector3dVec reconstVec;  // reconstruction of every pixel by cam1
    Vector3dVec reconstRotVec;  // reconstVec rotated into the second frame
    
    Eigen::Vector2d epipole;  // projection of the first camera center onto the second camera
    Vector2dVec pinfVec;  // projection of reconstRotVec by cam2
    
    // discretized version
    //TODO remove
    cv::Point epipolePx;  
    std::vector<cv::Point> pinfPxVec;
    
    std::vector<Polynomial2> epipolarVec;  // the epipolar curves represented by polynomial functions
    
    Mat8u errorBuffer;
    cv::Mat_<int> tableauLeft, tableauRight;
    cv::Mat_<int> tableauTop, tableauBottom;
    Mat8u smallDisparity;
    
    StereoParameters params;
};

