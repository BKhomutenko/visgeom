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
   
MonoOdometry::MonoOdometry(const ptree & params):
    _xiBaseCam( readTransform(params.get_child("xi_base_camera")) ),
    _sgmParams(params.get_child("stereo_parameters")),
    _camera( new EnhancedCamera(readVector<double>(params.get_child("camera_params")).data()) ),
    sparseOdom(_camera, _xiBaseCam),
    motionStereo(_camera, _camera, params.get_child("stereo_parameters")),
    localizer(5, _camera),
    state(STATE_BEGIN)
{
    localizer.setVerbosity(0);
    localizer.setXiBaseCam(_xiBaseCam);
}
    
MonoOdometry::~MonoOdometry() 
{
    delete _camera;
}
   
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
        _xiGlobal = _xiLocal = Transf(0, 0, 0, 0, 0, 0);
        imageVec.emplace_back(imageNew.clone());
        transfVec.push_back(_xiLocal);
        motionStereo.setBaseImage(imageNew);
        sparseOdom.feedData(imageNew, _xiLocal);
        localizer.setBaseImage(imageNew);
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
            EnhancedSgm sgm(getCameraMotion(), _camera, _camera, _sgmParams);
            //compute Sgm stereo
            sgm.computeStereo(imageVec.back(), imageNew, depth);
            
            //go to the next state
            localizer.setDepth(depth);
            state = STATE_READY;
            
            //Check the localization eigenvalues
            
            localizer.setTargetImage(imageNew);
            double epsilon = 0.00;
            Transf xiPhoto = localizer.computePose(_xiLocal).compose(Transf(epsilon, -epsilon, epsilon, epsilon, -epsilon, -epsilon));
            Transf xiMI = localizer.computePoseMI(_xiLocal);
            cout << "xiPhoto" << endl << xiPhoto << endl;
            cout << "xiMI" << endl << xiMI << endl;
            cout << "xiFeature" << endl << _xiLocal << endl;
            array<double, 6>  eigen1 = localizer.covarianceEigenValues(1, _xiLocal, false);
            for (int i = 0; i < 6; i++)
            {
                cout << setw(16) << eigen1[i];
            }
            cout << endl;
        }
        break; 
    case STATE_READY:
        
        localizer.setTargetImage(imageNew);
        _xiLocal = localizer.computePose(_xiLocal);
//        array<double, 6>  eigen1 = localizer.covarianceEigenValues(1, _xiLocal, false);
//            for (int i = 0; i < 6; i++)
//            {
//                cout << setw(16) << eigen1[i];
//            }
//            cout << endl;
        cout << _xiGlobal.compose(_xiLocal) << endl << endl;
        
        depth = motionStereo.compute(getCameraMotion(), imageNew, depth);
        depth.filterNoise();
        
        localizer.setDepth(depth);
        //TODO check whether a new keyframe is needed
        
        if (isNewKeyframeNeeded())
        {
            pushKeyFrame(imageNew);
        }
        
        break;
    }
}

void MonoOdometry::pushKeyFrame(const Mat8u & imageNew)
{
    _xiGlobal = _xiGlobal.compose(_xiLocal);
    transfVec.push_back(_xiLocal);
    _xiLocal = Transf(0, 0, 0, 0, 0, 0);
    imageVec.push_back(imageNew.clone());
    localizer.setBaseImage(imageNew);
    //Project the depth forward
    depth = depth.wrapDepth(_xiLocal);
}

bool MonoOdometry::isNewKeyframeNeeded()
{
    return (_xiLocal.trans().norm() > MAX_DIST or _xiLocal.rot().norm() > MAX_ANGLE);
}

Transf MonoOdometry::getCameraMotion() const
{
    return _xiBaseCam.inverseCompose(_xiLocal).compose(_xiBaseCam);
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
    

