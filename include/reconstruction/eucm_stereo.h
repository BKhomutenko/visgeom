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
Semi-global block matching algorithm for non-rectified images
NOTE:
(u, v) is an image point 
(x, y) is a depth map point 
*/

#pragma once

#include "std.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "camera/eucm.h"

#include "utils/scale_parameters.h"
#include "reconstruction/stereo_misc.h"

struct StereoParameters : public ScaleParameters
{
    int dispMax = 48;
  
    int maxError = 25;
    int maxBias = 10;
    
    int verbosity = 0;

    //multi-hypoteses
    int hypMax = 1;
    int maxHypDiff = 10;
    
    //for descriptor matching
    int flawCost = 7;
};

//TODO separate the triangulation from this class
//TODO change EnhancedCamera to ICamera

class EnhancedStereo
{
public:
    
    EnhancedStereo(const EnhancedCamera * cam1, const EnhancedCamera * cam2,
            const StereoParameters & parameters) :
            // initialize members
            params(parameters),
            camera1(cam1->clone()),
            camera2(cam2->clone()),
            Transform12(1, 0, 0, 0, 0, 0)
    { 
    }
    
    ~EnhancedStereo()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
    }
    
    void setTransformation(const Transformation<double> & T12) { Transform12 = T12; }
    
    bool triangulate(double u1, double v1, double u2, double v2, Vector3d & X) const;
    
    // returns the distance between corresponding camera and the point
    double triangulate(double u1, double v1, double u2, double v2, CameraIdx camIdx = CAMERA_1) const;
    
    vector<int> compareDescriptor(const vector<uint8_t> & desc, const vector<uint8_t> & sampleVec) const;
    
    //FIXME potentially put into misc file
    vector<int> findLocalMinima(const vector<int> & dataVec) const;
protected:
    StereoParameters params;
    
    Transformation<double> Transform12;  // pose of camera 2 wrt camera 1
    EnhancedCamera *camera1, *camera2;
};

