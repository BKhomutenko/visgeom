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
#include "localization/sparse_odom.h"

struct MappingParameters
{
    double distThreshSq = 0.25;
    double angularThreshSq = 0.25;
};


struct Frame
{
    Mat8u img;
    int parent;  //so far map is a tree structure
    Transf xi; //the position is defined in the global frame
};

//Speed estimation and extrapolation are to be added
//Assumed that the data arrives in the chronological order


class PhotometricMapping
{
public:
    PhotometricMapping(const MappingParameters & params);
    
    void feedOdometry(const Transf & xi);
    
    void feedImage(const Mat8u & img);
    
    void reInit(const Transf & xi);
    
    int select(); //-1 means that there is no matching frame
    
    void newFrame(const Mat8u & img);
    Transf localizeMI(const Mat8u & img) const;
    Transf localizePhoto(const Mat8u & img) const;
    
    Transf getCameraMotion() const;
    
private:
    enum State {ST_INIT, ST_SELECT, ST_TMP_FRAME, ST_MAP_FRAME};
    MappingParameters _params;
    State _state;
    
    int _activeFrameIdx;
    Frame _tmpFrame;
    vector<Frame> _frameVec;

    Transf _xi; //current estimation
    Transf _xiOdom; //the last odometry measure
    Transf _xiOdomImage; //the last odometry before the last image
    
    Transf _xiBaseCamera;
    
    DepthMap depth;
};







