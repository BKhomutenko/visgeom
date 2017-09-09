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
Relative camera pose estimation based on photometric error and depth map
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"
#include "ceres.h" 

#include "geometry/geometry.h"
#include "projection/generic_camera.h"

//TODO make a parameter structure
//const double MIN_INIT_DIST = 0.25;   // minimal distance traveled befor VO is used
//const double MIN_STEREO_BASE = 0.25; // minimal acceptable stereo base
 // Feature matching threshold

class SparseOdometry
{
public:
    SparseOdometry(const ICamera * cameraPtr,
                Transf xiBaseCam, double minDist = 0.):
    xiBaseCam(xiBaseCam),
    depthState(STATE_EMPTY),
    keyframeState(STATE_EMPTY),
    camera(cameraPtr->clone()),
//    detector(25, 3, 1),
    numRansacPoints(2),
    MIN_STEREO_BASE(minDist),
    _g(0)
    { }
        
    ~SparseOdometry() { delete camera; }
    
    void feedData(const Mat8u & imageNew, const Transf xiOdomNew);
    
      
    double computeTransfSparse(const Vector3dVec & xVec1, const Vector3dVec & xVec2, 
            const Vector2dVec & pVec2, const vector<double> & sizeVec, const Transf xiOdom, Transf & xiOut, bool report = false);
    
    void ransacNPoints(const Vector3dVec & cloud1,
        const Vector3dVec & cloud2, const Vector2dVec & ptVec2,
        const vector<double> & sizeVec,
        const Transf xiOdom, vector<bool> & inlierMask);
    
    void motionMatchesFilter(const Vector3dVec & cloud1,
        const Vector3dVec & cloud2, const Vector2dVec & ptVec2,
        const Transf xiOdom, vector<bool> & inlierMask);
    
    Vector2dVec reprojectPoints(const Vector3dVec & cloud1,
    const Vector3dVec & cloud2, const Transf dxi);
    
    const Transf & getIncrement() const { return xiIncr; }
    const Transf & getIntegrated() const { return xiLocal; }
    
private:
    
    //TODO TEST
    Vector2dVec keypointVec1;
    Vector2dVec keypointVec2;
    
//    vector<KeyPoint> keypointVec1;
//    vector<KeyPoint> keypointVec2;
    
    Mat32f desc1, desc2;
    
    Mat8u imageOld;
    
    ICamera * camera;
    
    // 3 points fully constraint the transformation
    // 1 or 2 points rely on the odometry estimation
    const int numRansacPoints; 
    
    
    // state
    Transf xiLocal; // integrated VO path
    Transf xiIncr; // the last pose increment
    Transf xiOdom; // the last WO measurement
    const Transf xiBaseCam; // extrinsic calibration
    Transf xi5, xiTmp;

    enum DataState {STATE_EMPTY, STATE_READY};

    DataState depthState;
    DataState keyframeState;
    BRISK detector;
//    cv::ORB detector;
    mt19937 _g;
    
    const int distThresh = 2500; //for descriptor comparison
    const double MIN_STEREO_BASE;  // minimal acceptable stereo base
};



            
