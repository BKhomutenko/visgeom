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
#include "curve_rasterizer.h"
#include "depth_map.h"
#include "eucm_epipolar.h"

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
    int maxBias = 10;
    
    int verbosity = 0;
    int maxDepth = 100;
    // precomputed parameters
    int u0, v0;
    int uMax, vMax;
    int xMax, yMax;
    int halfBlockSize;
    
    // to be called before using
    void init()
    {
        u0 = uMargin + scale; 
        v0 = vMargin + scale;
        
        if (width > 0) uMax = u0 + width;
        else uMax = imageWidth - uMargin - scale;
        
        if (height > 0) vMax = v0 + height;
        else vMax = imageHeight - vMargin - scale;
        
        xMax = x(uMax) + 1;
        yMax = y(vMax) + 1;
        
        halfBlockSize =  scale / 2; 
    }
    
    // from image to small disparity coordiante transform
    int x(double u) { return round((u - u0) / scale); }
    int y(double v) { return round((v - v0) / scale); }
    
    // from small disparity to image coordiante transform    
    int u(int x) { return x * scale + u0; }
    int v(int y) { return y * scale + v0; }
    
    // from small disparity to image block corner transform 
//    int uCorner(int u) { return u * scale + u0; }
//    int vCorner(int v) { return v * scale + v0; }
};

class EnhancedStereo
{
public:
    enum CameraIdx {CAMERA_1, CAMERA_2};
    
    EnhancedStereo(Transformation<double> T12, const EnhancedCamera * cam1,
            const EnhancedCamera * cam2, const StereoParameters & stereoParams) :
            // initialize members
            Transform12(T12), 
            camera1(cam1->clone()),
            camera2(cam2->clone()),
            params(stereoParams),
            epipolar(T12, cam1, cam2, 2500)
    { 
        assert(params.dispMax % 2 == 0);
        params.init();
        init(); 
    }
    
    ~EnhancedStereo()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
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
    }
    
    // An interface function
    void comuteStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depthMap);
    
    // An interface function
    void comuteStereo(const Mat8u & img1, const Mat8u & img2, Mat32f & depthMat);
    
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
    void computeEpipolarIndices();
    
    // TODO remove from this class
    // draws an epipolar line  on the right image that corresponds to (x, y) on the left image
    void traceEpipolarLine(int u, int v, Mat8u & out, CameraIdx camIdx);
    
    //// DYNAMIC PROGRAMMING
    void createBuffer();
       
    // fill up the error buffer using 2*S-1 pixs along epipolar lines as local desctiprtors
    void computeCurveCost(const Mat8u & img1, const Mat8u & img2);
    
    void computeDynamicProgramming();
    
    void computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost);
    
    void reconstructDisparity();  // using the result of the dynamic programming
    
    // TODO implement
    void upsampleDisparity(const Mat8u & img1, Mat8u & disparityMat);
    
    //// MISCELLANEOUS
    
    // index of an object in a linear array corresponding to pixel [row, col] 
    int getLinearIndex(int x, int y) { return params.xMax*y + x; }
      
    CurveRasterizer<int, Polynomial2> getCurveRasteriser1(int idx);
    CurveRasterizer<int, Polynomial2> getCurveRasteriser2(int idx);
    
    // reconstruction
    bool triangulate(double u1, double v1, double u2, double v2, Vector3d & X);
    void computeDepth(Mat32f & distanceMat);
    
    //TODO put generatePlane elsewhere
    void generatePlane(Transformation<double> TcameraPlane, 
            Mat32f & distance, const Vector3dVec & polygonVec);
            
    //TODO put generatePlane elsewhere        
    void generatePlane(Transformation<double> TcameraPlane, 
            DepthMap & distance, const Vector3dVec & polygonVec);
            
    double computeDepth(int x, int y);
    bool computeDepth(int x, int y, double & dist, double & sigma);
    
private:
    EnhancedEpipolar epipolar;
    
    Transformation<double> Transform12;  // pose of camera 2 wrt camera 1
    EnhancedCamera *camera1, *camera2;
   
    std::vector<bool> maskVec;
    
    Vector2dVec pointVec1;  // the depth points on the image 1
    Vector3dVec reconstVec;  // reconstruction of every pixel by cam1
    Vector3dVec reconstRotVec;  // reconstVec rotated into the second frame
    
    bool epipoleInverted1, epipoleInverted2;
    Vector2d epipole1, epipole2;  // projection of the first camera center onto the second camera
    Vector2dVec pinfVec;  // projection of reconstRotVec by cam2
    
    // discretized version
    Vector2iVec pointPxVec1;
    Vector2i epipolePx1, epipolePx2;  
    Vector2iVec pinfPxVec;
    
    Mat8u errorBuffer;
    Mat32s tableauLeft, tableauRight; //FIXME check the type through the code
    Mat32s tableauTop, tableauBottom;
    Mat8u smallDisparity;
    
    StereoParameters params;
};

