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
SLAM system
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"

#include "reconstruction/depth_map.h"
#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "projection/eucm.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/eucm_sgm.h"
#include "localization/sparse_odom.h"
#include "localization/photometric.h"

struct MappingParameters
{
    MappingParameters(const ptree & params) 
    {
        for (auto & item : params)
        {
            const string & pname = item.first;
            if (pname == "new_kf_distance_thresh") distThreshSq = pow(item.second.get_value<double>(), 2);
            else if (pname == "new_kf_angular_thresh") angularThreshSq = pow(item.second.get_value<double>(), 2);
            else if (pname == "min_init_dist")  minInitDist = item.second.get_value<double>();
            else if (pname == "min_stereo_base") minStereoBase = item.second.get_value<double>();
            else if (pname == "dist_thresh") maxDistance = item.second.get_value<double>();
            else if (pname == "normalize_scale") normalizeScale = item.second.get_value<bool>();
        }
        
    }
    
    MappingParameters() {}
    //defines the conditions for the new frame instantiation
    double distThreshSq = 0.03;
    double angularThreshSq = 0.25;
    
    //minimum distance to initialize the depth map and the whole system
    double minInitDist = 0.1;
    
    //if the stereo base is below the stereo is not computed
    double minStereoBase = 0.07;
    
    //beyond this distance the points are not used for the localization
    double maxDistance = 5;
    bool normalizeScale = true;
};


struct Frame
{
    Mat8u img;
    Transf xi; //the position is defined in the global frame
};

//Speed estimation and extrapolation are to be added
//Assumed that the data arrives in the chronological order


class PhotometricMapping
{
public:
    PhotometricMapping(const ptree & params);
    
    //true if the frame has been inserted
    bool constructMap(const Transf & xi, const Mat8u & img);
    
    void feedOdometry(const Transf & xi);
    
    void feedImage(const Mat8u & img);
    
    void pushInterFrame(const Mat8u & img);
    
    void reInit(const Transf & xi);
    
    int selectMapFrame(const Transf & xi, const double K = 4); //-1 means that there is no matching frame
    
    Transf localizeMI(); //localizes the inter frame wrt _frameVec[_mapIdx]
    Transf localizePhoto(const Mat8u & img); //localizes the image wrt interFrame
    
    Transf getCameraMotion(const Transf & xi) const;
    
    bool checkDistance(const Transf & xi, const double K = 1.) const;
    bool checkDistance(const Transf & xi1, const Transf & xi2, const double K = 1.) const;
    
    void improveStereo(const Mat8u & img);
//private:

    enum State {MAP_BEGIN, MAP_INIT, MAP_LOCALIZE, MAP_SLAM};
    enum WarningType {WARNING_NONE, WARNING_SCALE, WARNING_ROTATION};
    
    
    
    MappingParameters _params;
    
    EnhancedCamera * _camera;
    //state variables
    State _state;
    WarningType _warningState;
    
    int _mapIdx; //currently active map frame
    bool _odomInit;
    Frame _interFrame;
    vector<Frame> _frameVec;
    DepthMap _depth;
    Transf _xiLocal; //current base pose estimation in the local frame
    Transf _xiLocalOld; //for VO scale rectification
    Transf _xiOdom; //the last odometry measure
    Transf _zetaOdom; //the last odometry measure
    Transf _xiOdomImage; //the last odometry before the last image
    
    Transf _xiBaseCam;
    
    //utils
    //are used to create an Sgm object to init the keyframe
    SgmParameters _sgmParams; 
    
    //used to initialize the first transformation
    SparseOdometry _sparseOdom;
    
    //gradually improves the depth map
    MotionStereo _motionStereo;
    
    ScalePhotometric _localizer;
};







