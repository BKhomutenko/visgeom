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
const double MIN_INIT_DIST = 0.25;   // minimal distance traveled befor VO is used
const double MIN_STEREO_BASE = 0.15; // minimal acceptable stereo base
const int distThresh = 100;

class SparseOdometry
{
public:
    SparseOdometry(const ICamera * cameraPtr,
                Transf xiBaseCam):
    xiBaseCam(xiBaseCam),
    depthState(STATE_EMPTY),
    keyframeState(STATE_EMPTY),
    camera(cameraPtr->clone()),
    detector(40, 0, 2) 
    
    { }
        
    ~SparseOdometry() { delete camera; }
    
    void feedData(const Mat8u & imageNew, const Transf xiOdomNew);
    
      
    double computeTransfSparse(const Vector3dVec & xVec1, const Vector3dVec & xVec2, 
            const Vector2dVec & pVec2, const Transf xiOdom, Transf & xiOut, bool report = false);
    
    void ransacFivePoint(const Vector3dVec & cloud1,
        const Vector3dVec & cloud2, const Vector2dVec & ptVec2,
        const Transf xiOdom, vector<bool> & inlierMask);
        
private:
    vector<KeyPoint> keypointVec1;
    vector<KeyPoint> keypointVec2;
    
    Mat32f desc1, desc2;
    
    ICamera * camera;
    
    // state
    Transf xiLocal; // pose wrt actual keyframe
    Transf xiOdom; // the last WO measurement
    const Transf xiBaseCam; // extrinsic calibration
    Transf xi5, xiTmp;

    enum DataState {STATE_EMPTY, STATE_READY};

    DataState depthState;
    DataState keyframeState;
    BRISK detector;
//    cv::ORB detector;
    
};



            
