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

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/eucm_sgm.h"
#include "localization/photometric.h"
#include "localization/sparse_odom.h"

//TODO make a parameter structure
const double MIN_INIT_DIST = 0.25;   // minimal distance traveled befor VO is used
const double MIN_STEREO_BASE = 0.05; // minimal acceptable stereo base

//constants to create new keyframes, TODO put elsewhere
const double MAX_DIST = 0.5;
const double MAX_ANGLE = M_PI / 4;

//TODO parameters initialization via .json
class MonoOdometry
{
public:
    MonoOdometry(const ptree & params);
        
    ~MonoOdometry();
    
    void feedWheelOdometry(const Transf xiOdomNew);
    void feedImage(const Mat8u & imageNew);
    
    
    Transf getCameraMotion() const;
    
//private:

    void pushKeyFrame(const Mat8u & imageNew);
    bool isNewKeyframeNeeded();
    
    EnhancedCamera * _camera;
    // memory
    vector<Mat8u> imageVec; //.back() is the actual key frame
    vector<Transf> transfVec;
    vector<Matrix6d> poseCovarVec;
    
    // state
    // pose wrt current keyframe
    Transf _xiLocal; 
    Transf _xiGlobal;
    // the last WO measuremen to compute odometry increment
    Transf _xiOdom; 
    
    // camera position wrt wheel odometry frame
    Transf _xiBaseCam; 
    
    // related to the key frame
    DepthMap depth; 

    enum DataState {STATE_BEGIN, STATE_SPARSE_INIT, STATE_READY};

    DataState state;
    
    //utils
    //are used to create an Sgm object to init the keyframe
    SgmParameters _sgmParams; 
    
    //used to initialize the first transformation
    SparseOdometry sparseOdom;
    
    //gradually improves the depth map
    MotionStereo motionStereo;
    
    ScalePhotometric localizer;
};



