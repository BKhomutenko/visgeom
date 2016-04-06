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
Stereo vision definition.
*/

#pragma once
//STL
#include <vector>
#include <array>
#include <algorithm>

//Eigen
#include <Eigen/Eigen>

//OpenCV
#include <opencv2/opencv.hpp>

#include "geometry.h"
#include "eucm.h"
#include "curve_rasterizer.h"

using namespace std;

typedef CurveRasterizer<Polynomial2> EpipolarRasterizer;

class StereoEUCM
{
public:
    StereoEUCM(Transformation<double> T12, int imageWidth, int imageHeight,
            const double * params1, const double * params2)
            : Transform12(T12), 
            cam1(imageWidth, imageHeight, params1), 
            cam2(imageWidth, imageHeight, params2),
            dispMax(64), verticalMargin(100), blockSize(5),
            lambdaStep(6), lambdaJump(20)
            { init(); }

    void setTransformation(Transformation<double> T12) 
    { 
        Transform12 = T12;
        initAfterTransformation();
    } 
    
    void init()
    {
        computeReconstructed();
        initAfterTransformation();
    }
    
    void initAfterTransformation()
    {
        computeEpipole();
        computeRotated();
        computePinf();
        computeEpipolarCurves();
    }
    
    // **EPIPOLAR GEOMETRY**
    
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
    void traceEpipolarLine(cv::Point pt, cv::Mat_<uint8_t> & out);
    
    
    // **DYNAMIC PROGRAMMING**
    // fills up the output with photometric errors between the val = I1(pi) and 
    // the values from I2 on the epipolar curve
    void comuteStereo(const cv::Mat_<uint8_t> & img1, 
            const cv::Mat_<uint8_t> & img2,
            cv::Mat_<uint8_t> & disparity);
    
    void computeCost(const cv::Mat_<uint8_t> & img1, const cv::Mat_<uint8_t> & img2);
    
    void computeDynamicProgramming();
    
    void computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost);
    
    void reconstructDisparity();  // using the result of the dynamic programming
    
    void upsampleDisparity(const cv::Mat_<uint8_t> & img1, cv::Mat_<uint8_t> & disparity);
    
    void getLinearIdx(int row, int col) { return cam1.width*row + col; }
    
    void uSmall(int u) { return (u - u0) / blockSize; }
    
    void uBig(int u) { return u * blockSize + blockSize/2 + u0; }
    
    void vSmall(int v) { return (v - v0) / blockSize; }
    
    void vBig(int v) { return v * blockSize + blockSize/2 + v0; }
    
private:
    Transformation<double> Transform12;  // pose of the first to the second camera
    EnhancedCamera cam1, cam2;
   
    vector<Eigen::Vector3d> reconstVec;  // reconstruction of every pixel by cam1
    vector<Eigen::Vector3d> reconstRotVec;  // reconstVec rotated into the second frame
    
    Eigen::Vector2d epipole;  // projection of the first camera center onto the second camera
    vector<Eigen::Vector2d> pinfVec;  // projection of reconstRotVec by cam2
    
    // discretized version
    cv::Point epipolePx;  
    vector<cv::Point> pinfPxVec;
    
    vector<Polynomial2> epipolarVec;  // the epipolar curves represented by polynomial functions
    
    int dispMax; // maximum shift along the epipolar line
    int u0, v0, uMax, vMax;  // the active rectangle 
    int blockSize;
    
    int lambdaStep, lambdaJump;  // costs for disparity changes
    
    cv::Mat_<uint8_t> errorBuffer;
    cv::Mat_<int> tableauLeft, tableauRight;
    cv::Mat_<int> tableauTop, tableauBottom;
    cv::Mat_<uint8_t> smallDisparity;
};

