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
#include "camera/generic_camera.h"
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
                Transformation<double> xiBaseCam,
                MotionStereoParameters params):
    xiBaseCam(xiBaseCam),
    motionStereo(camera, camera, params),
    photometricLocalizer(5, camera),
    state(EMPTY),
    woState(EMPTY)
    {
        photometricLocalizer.setVerbosity(0);
    }
        
    ~MonoOdometry() {}
    
    void feedImage(const Mat8u & imageNew, const double t);
    
    void feedWheelOdometry(const Transformation<double> xiOdomNew, const double t);
    
private:

    void initFirstKeyFrame(const Mat8u & imageNew, const double t);
    
    //xiLocal is supposed to be up to date
    void createKeyFrame(const Mat8u & imageNew, const double t);
    
    void computeVisualOdometry(const Mat8u & imageNew, const double t);
    
    void refineDepth(const Mat8u & imageNew);
    
    enum OdomState {EMPTY,                 // no data
                    INIT_FIRST_FRAME,      // no good depth estimation
                    ACTIVE,                 
                    LOST};
    // memory
    vector<Mat8u> imageVec;
    vector<Transformation<double>> transfoVec;
    vector<Matrix6d> poseCovarVec;
    
    // state
    Transformation<double> xiLocal; // pose wrt actual keyframe
    Transformation<double> xiOdom; // the last WO measurement
    Transformation<double> xiBaseCam; // extrinsic calibration
    DepthMap depthMap; 
    OdomState state, woState;
    double tLocal, tOdom;
    
    //utils
    MotionStereo motionStereo;
    ScalePhotometric photometricLocalizer;
    //TODO add WO sigma
    //TODO add a threshold for new the KFr instantiation 
};



