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
    cout << "verbosity " << _sgmParams.verbosity << endl; 
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

void MonoOdometry::setDepth(const Mat32f & imageDepth)
{
    if (depth.empty())
    {
        depth = DepthMap(_camera, _sgmParams);
    }
    for (int y = 0; y < depth.yMax; y++)
    {
        for (int x = 0; x < depth.xMax; x++)
        {
            int u = _sgmParams.uConv(x);
            int v = _sgmParams.vConv(y);
            depth.at(x, y) = imageDepth(v, u);
            depth.sigma(x, y) = 0.1;
        }
    }
}

void MonoOdometry::feedImage(const Mat8u & imageNew)
{
    switch (state)
    {
    case STATE_BEGIN:
        //the starting point, position 0
        //TODO refactor
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
//        motionStereo.setBaseImage(imageNew);
//        depth = motionStereo.compute(getCameraMotion().inverse(), imageVec.back());
//        pushKeyFrame(imageNew);
        if (_xiLocal.trans().norm() >= MIN_INIT_DIST)
        {
            //compute accurate motion estimation
            sparseOdom.feedData(imageNew, _xiLocal);
            _xiLocal = sparseOdom.getIncrement();
            _xiLocalOld = _xiLocal;
            //create a keyframe
            pushKeyFrame(imageNew);
            
            //change the state
            state = STATE_READY;   
        }
        break; 
    case STATE_READY:
        localizer.setDepth(depth);
        localizer.setTargetImage(imageNew);
        
//        _xiLocal = localizer.computePose( _xiLocal.compose(Transf(0.01, 0.01, 0.01, 0.01, 0.01, 0.01)) );
        
        //scale factor normalization
        //TODO include the distance constraint in the optimization loop
        Transf zetaPrior = _xiLocalOld.inverseCompose(_xiLocal);
        _xiLocal = localizer.computePose( _xiLocal );
        Transf zeta = _xiLocalOld.inverseCompose(_xiLocal);
        double lengthPrior = zetaPrior.trans().norm();
        double length = zeta.trans().norm();
        double K = lengthPrior / length;
        if (length > 1e-2) //TODO put a meaningful threshold
        {
            
            if (K > 1.2) //too much divergence between WO and VO
            {
                cout << "WARNING : discard VO result, scale error" << endl;
                _xiLocal = _xiLocalOld.compose(zetaPrior);
                break;
            }
            zeta.trans() *= K;
            _xiLocal = _xiLocalOld.compose(zeta);
        }
        
        cout << "MOTION ESTIMATION :" << endl;
        cout << "    " << _xiLocal << endl;
        cout << "    K : " << K << endl;
        _xiLocalOld = _xiLocal;
//        array<double, 6>  eigen1 = localizer.covarianceEigenValues(1, _xiLocal, false);
//            for (int i = 0; i < 6; i++)
//            {
//                cout << setw(16) << eigen1[i];
//            }
//            cout << endl;
//        cout << _xiGlobal.compose(_xiLocal) << endl << endl;
        
        
        
        
        //TODO check whether a new keyframe is needed
        
        if (isNewKeyframeNeeded())
        {
            pushKeyFrame(imageNew);
        }
        else
        {
            depth = motionStereo.compute(getCameraMotion(), imageNew, depth);
            depth.filterNoise();
        }
        
        break;
    }
   
    cout << "composed :" << endl;
    cout << "    " << _xiGlobal.compose(_xiLocal) << endl;
}

void MonoOdometry::pushKeyFrame(const Mat8u & imageNew)
{

    //set Sgm stereo base
    EnhancedSgm sgm(getCameraMotion().inverse(), _camera, _camera, _sgmParams);
    //compute Sgm stereo 1-0
    DepthMap depthNew;
    sgm.computeStereo(imageNew, imageVec.back(), depthNew);
    
    
//    //set Sgm stereo base
//    EnhancedSgm sgm(getCameraMotion(), _camera, _camera, _sgmParams);
//    //compute Sgm stereo 1-0
//    DepthMap depthNew;
//    sgm.computeStereo(imageVec.back(), imageNew, depthNew);
//    
    
    
    if (state == STATE_READY)
    {
        // project the prior depth estimation       
        depth = depth.wrapDepth(_xiLocal);
        depth.merge(depthNew);
//        depth = depthNew;
    }
    else
    {
        depth = depthNew;
    }
    
    //store the motion estimation
    _xiGlobal = _xiGlobal.compose(_xiLocal);
    transfVec.push_back(_xiLocal);
    
    //reset the local position
    _xiLocal = Transf(0, 0, 0, 0, 0, 0);
    imageVec.emplace_back(imageNew.clone());
    localizer.setBaseImage(imageNew);
    motionStereo.setBaseImage(imageNew);
    //Project the depth forward
    
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
    

