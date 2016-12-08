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

    
void MonoOdometry::feedData(const Mat8u & imageNew, const Transf xiOdomNew)
{
    switch (keyframeState)
    {
    case STATE_EMPTY:
        _xiOdom = xiOdomNew;
        motionStereo.setBaseImage(imageNew);
        keyframeState = STATE_READY;   
        break;     
    case STATE_READY:
        //integrate WO
        _xiLocal = _xiLocal.composeInverse(_xiOdom).compose(xiOdomNew);
        
        if (depthState == STATE_READY)
        {
            computeVisualOdometry(imageNew);
            //TODO check conditioning
        }
        
        refineDepth(imageNew);
    }
}
    
//void MonoOdometry::initFirstKeyFrame(const Mat8u & imageNew)
//{
//    imageVec.emplace_back();
//    imageNew.copyTo(imageVec.back());
//    transfoVec.emplace_back();
//    motionStereo.setBaseImage(imageNew);
//    
//    //update the state
//    xiLocal = Transf();
//    tLocal = t;
//    
//    //TODO rename to setBaseImage
//	photometricLocalizer.computeBaseScaleSpace(imageNew);
//	
//    //TODO treat covariance
//    //poseCovarVec.emplaceBack();
//}
    
//void MonoOdometry::createKeyFrame(const Mat8u & imageNew)
//{
//    // project depth forward
//    
//    // insert new instances into the map
//    
//    // update the states
//}

void MonoOdometry::computeVisualOdometry(const Mat8u & imageNew)
{
    Transf xiCam12 = _xiBaseCam.inverseCompose(_xiLocal.compose(_xiBaseCam));;
    
    photometricLocalizer.depth() = depthMap; //TODO make a method setDepth()
    
    photometricLocalizer.computePose(imageNew, xiCam12); //TODO optimize directly xiLocal
                            //TODO get the covariance matrix
                            //TODO introduce the odometry prior with its propoer covariance
    
    _xiLocal = _xiBaseCam.compose(xiCam12.composeInverse(_xiBaseCam));
}

void MonoOdometry::refineDepth(const Mat8u & imageNew)
{
    Transf xiCam12 = _xiBaseCam.inverseCompose(_xiLocal).compose(_xiBaseCam);
    double base = xiCam12.trans().norm();
    if (base < MIN_STEREO_BASE) return; 
    motionStereo.validateDepth(xiCam12, imageNew, depthMap);
    if (depthState == STATE_EMPTY and base > MIN_INIT_DIST) depthState = STATE_READY;
}
    

