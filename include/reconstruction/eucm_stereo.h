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
            dispMax(64), verticalMargin(100), blockSize(5)
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
    
    // fills up the output with photometric errors between the val = I1(pi) and 
    // the values from I2 on the epipolar curve
    void computeCost(const cv::Mat_<uint8_t> & img1, const cv::Mat_<uint8_t> & img2, cv::Mat_<uint8_t> & out);
    
    void computeDynamicProgramming(const cv::Mat_<uint8_t> & costMat, cv::Mat_<int> & out);
    
    // to visualize the epipolar lines
    void traceEpipolarLine(cv::Point pt, cv::Mat_<uint8_t> & out);
    
private:
    Transformation<double> Transform12;  // pose of the first to the second camera
    EnhancedCamera cam1, cam2;
   
    vector<Eigen::Vector3d> reconstVec;  // reconstruction of every pixel by cam1
    vector<Eigen::Vector3d> reconstRotVec;  // reconstVec rotated into the second frame
    
    Eigen::Vector2d epipole;
    vector<Eigen::Vector2d> pinfVec;  // projection of reconstRotVec by cam2
    
    cv::Point epipolePx;  // Px means that it is discrete
    vector<cv::Point> pinfPxVec;
    
    vector<Polynomial2> epipolarVec;  // the epipolar curves represented by polynomial functions
    
    int dispMax; //maximum shift along the epipolar
    int verticalMargin;  //number of rows that we skip from the top and bottom of the image
    int blockSize;
    
};

