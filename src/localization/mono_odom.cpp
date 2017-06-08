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
#include "projection/generic_camera.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"

#include "localization/sparse_odom.h"
   
void MonoOdometry::feedWheelOdometry(const Transf xiOdomNew)
{
    //integrate
    _xiLocal = _xiLocal.composeInverse(_xiOdom).compose(xiOdomNew);
    
    //refresh
    _xiOdom = xiOdomNew;
}

void MonoOdometry::feedImage(const Mat8u & imageNew)
{
    switch (state)
    {
    case STATE_BEGIN:
        //the starting point
        _xiLocal = Transf(0, 0, 0, 0, 0, 0);
        imageVec.emplace_back(imageNew.clone());
        motionStereo.setBaseImage(imageNew);
        sparseOdom.feedData(imageNew, _xiLocal);
        state = STATE_SPARSE_INIT;   
        break;    
    case STATE_SPARSE_INIT:
        //wait until the distance is sucfficient to do the saprce initialization
        if (_xiLocal.trans().norm() >= MIN_INIT_DIST)
        {
            //compute accurate motion estimation
            sparseOdom.feedData(imageNew, _xiLocal);
            _xiLocal = sparseOdom.getIncrement();
            //compute the initial depth estimation
            
            //set Sgm stereo base
            Transf stereoBase = _xiBaseCam.inverseCompose(_xiLocal).compose(_xiBaseCam);
            EnhancedSgm sgm(_xiLocal, _camera, _camera, _sgmParams);
            
            //compute Sgm stereo
            sgm.computeStereo(imageVec.back(), imageNew, depth);
            //go to the next state
            state = STATE_READY;
        }
        break; 
    case STATE_READY:
        //TODO add code here
        break;
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

//void MonoOdometry::computeVisualOdometry(const Mat8u & imageNew)
//{
//    Transf xiCam12 = _xiBaseCam.inverseCompose(_xiLocal.compose(_xiBaseCam));;
//    
//    photometricLocalizer.depth() = depthMap; //TODO make a method setDepth()
//    
//    photometricLocalizer.computePose(imageNew, xiCam12); //TODO optimize directly xiLocal
//                            //TODO get the covariance matrix
//                            //TODO introduce the odometry prior with its propoer covariance
//    
//    _xiLocal = _xiBaseCam.compose(xiCam12.composeInverse(_xiBaseCam));
//}

//void MonoOdometry::refineDepth(const Mat8u & imageNew)
//{
//    Transf xiCam12 = _xiBaseCam.inverseCompose(_xiLocal).compose(_xiBaseCam);
//    double base = xiCam12.trans().norm();
//    if (base < MIN_STEREO_BASE) return; 
//    depthMap = motionStereo.compute(xiCam12, imageNew, depthMap);
//    if (depthState == STATE_EMPTY and base > MIN_INIT_DIST) depthState = STATE_READY;
//}
    

