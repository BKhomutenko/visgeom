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
#include "localization/photometric.h"

//TODO make a parameter structure
const double MIN_INIT_DIST = 0.25;   // minimal distance traveled befor VO is used
const double MIN_STEREO_BASE = 0.05; // minimal acceptable stereo base

class MonoOdometry
{
public:
    MonoOdometry(const EnhancedCamera * camera,
                Transf xiBaseCam,
                MotionStereoParameters params):
    _xiBaseCam(xiBaseCam),
    sparseOdom(camera, xiBaseCam),
    motionStereo(camera, camera, params),
    photometricLocalizer(5, camera),
    depthState(STATE_EMPTY),
    keyframeState(STATE_EMPTY)
    {
        photometricLocalizer.setVerbosity(0);
    }
        
    ~MonoOdometry() {}
    
    void feedData(const Mat8u & imageNew, const Transf xiOdomNew);
    
//private:
    //xiLocal is supposed to be up to date
    void computeVisualOdometry(const Mat8u & imageNew);
    
    void refineDepth(const Mat8u & imageNew);

    // memory
    vector<Mat8u> imageVec; //.back() is the actual key frame
    vector<Transf> transfVec;
    vector<Matrix6d> poseCovarVec;
    
    // state
    Transf _xiLocal; // pose wrt actual keyframe
    Transf _xiOdom; // the last WO measurement
    Transf _xiBaseCam; // extrinsic calibration
    DepthMap depthMap; 

    enum DataState {STATE_BEGIN, STATE_SPARSE_INIT, STATE_READY};

    DataState depthState;
    DataState keyframeState;
    
    //utils
    SparseOdometry sparseOdom;
    MotionStereo motionStereo;
    ScalePhotometric photometricLocalizer;
};



