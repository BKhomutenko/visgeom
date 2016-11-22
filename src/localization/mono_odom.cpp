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

#include "localization/mono_odom.h"

#include "std.h"
#include "eigen.h"
#include "ocv.h"

#include "geometry/geometry.h"
#include "camera/generic_camera.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"

    
void MonoOdometry::feedImage(const Mat8u & imageNew, const double t)
{
    switch (state)
    {
    case EMPTY:
        initFirstKeyFrame(imageNew);
        state = INIT_FIRST_FRAME;
    case INIT_FIRST_FRAME:
        refineDepth(imageNew);
        if (xiLocal.trans().norm() > MIN_TRANS) state = ACTIVE;        
    case ACTIVE:
        if (false) // conditionning gets wors
        {
            createKeyFrame(imageNew);
        }
        else
        {
            computeVisualOdometry(imageNew, t);
            refineDepth(imageNew);
        }
    }
}

void MonoOdometry::feedWheelOdometry(const Transformation<double> xiOdomNew, const double t)
{
    Transformation<double> dxi = xiOdom.inverseCompose(xiOdomNew);
    
    //TODO treat the case when tLocal < tOdom ???
    
    if (tLocal > tOdom)
    {
        //interpolate
        double lambda = (t - tLocal) / (t - tOdom);
        dxi.scale(lambda);
    }
    
    //update state
    xiLocal = xiLocal.compose(dxi);
    xiOdom = xiOdomNew;
    tOdom = t;
    tLocal = t;
}
    
void MonoOdometry::initFirstKeyFrame(const Mat8u & imageNew)
{
    imageVec.emplace_back();
    imageNew.copyTo(imageVec.back());
    transfoVec.emplace_back();
    motionStereo.setBaseImage(imageNew);
    
    //TODO rename to setBaseImage
	photometricLocalizer.computeBaseScaleSpace(imageNew);
	
    //TODO treat covariance
    //poseCovarVec.emplaceBack();
}
    
void MonoOdometry::createKeyFrame(const Mat8u & imageNew)
{
    //TODO project depth forward
}

void MonoOdometry::computeVisualOdometry(const Mat8u & imageNew, const double t)
{
    Transformation<double> xiCam12;
    photometricLocalizer.depth() = depthMap; //TODO make a method setDepth()
    photometricLocalizer.computePose(imageNew, xiCam12);
    xiLocal = xiBaseCam.compose(xiCam12.composeInverse(xiBaseCam));
    tLocal = t;
}

void MonoOdometry::refineDepth(const Mat8u & imageNew)
{
    motionStereo.computeDepth(xiLocal, imageNew, depthMap);
}
    
//    enum OdomState {EMPTY, INIT_FIRST_FRAME, ACTIVE, LOST};
//    // memory
//    vector<Mat8u> imageVec;
//    vector<Transformation<double>> transfoVec;
//    vector<Matrix6d> poseCovarVec;
//    
//    // state
//    Transformation<double> xiLocal; // pose wrt actual keyframe
//    Transformation<double> xiOdom; // the last WO measurement
//    Transformation<double> xiBaseCam; // extrinsic calibration
//    DepthMap depthMap; 
//    OdomState state;
//    
//    //utils
//    MotionStereo motionStereo;
//    
//    //TODO add WO sigma
//    //TODO add a threshold for new the KFr instantiation 




