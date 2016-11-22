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
        initFirstKeyFrame(imageNew, t);
        state = INIT_FIRST_FRAME;
        break;
    case INIT_FIRST_FRAME:
        refineDepth(imageNew);
        if (xiLocal.trans().norm() > MIN_INIT_DIST) state = ACTIVE;    
        break;    
    case ACTIVE:
    
        // if some singular values of the covariance matrix after the localization is too low
        
        computeVisualOdometry(imageNew, t);
        
        if (false) // conditionning gets wors
        {
            createKeyFrame(imageNew, t);
        }
        else
        {
            refineDepth(imageNew);
        }
        break;
    }
}

void MonoOdometry::feedWheelOdometry(const Transformation<double> xiOdomNew, const double t)
{
    //TODO potentially dangerous code because two automata work simultaneously
    // verification's needed
    switch (woState)
    {
    case EMPTY:
        tOdom = t;
        tLocal = t;
        xiOdom = xiOdomNew; 
        woState = ACTIVE;   
        break;     
    case ACTIVE:
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
        break;
    }
}
    
void MonoOdometry::initFirstKeyFrame(const Mat8u & imageNew, const double t)
{
    imageVec.emplace_back();
    imageNew.copyTo(imageVec.back());
    transfoVec.emplace_back();
    motionStereo.setBaseImage(imageNew);
    
    //update the state
    xiLocal = Transformation<double>();
    tLocal = t;
    
    //TODO rename to setBaseImage
	photometricLocalizer.computeBaseScaleSpace(imageNew);
	
    //TODO treat covariance
    //poseCovarVec.emplaceBack();
}
    
void MonoOdometry::createKeyFrame(const Mat8u & imageNew, const double t)
{
    // project depth forward
    
    // insert new instances into the map
    
    // update the states
}

void MonoOdometry::computeVisualOdometry(const Mat8u & imageNew, const double t)
{
    Transformation<double> xiCam12 = xiBaseCam.inverseCompose(xiLocal.compose(xiBaseCam));;
    
    photometricLocalizer.depth() = depthMap; //TODO make a method setDepth()
    
    photometricLocalizer.computePose(imageNew, xiCam12); //TODO optimize directly xiLocal
                            //TODO get the covariance matrix
                            //TODO introduce the odometry prior with its propoer covariance
    
    xiLocal = xiBaseCam.compose(xiCam12.composeInverse(xiBaseCam));
    tLocal = t;
}

void MonoOdometry::refineDepth(const Mat8u & imageNew)
{
    if (xiLocal.trans().norm() < MIN_STEREO_BASE) return; 
    motionStereo.computeDepth(xiLocal, imageNew, depthMap);
}
    

