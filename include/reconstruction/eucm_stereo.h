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
#include "projection/eucm.h"

#include "reconstruction/scale_parameters.h"
#include "reconstruction/triangulator.h"
#include "reconstruction/epipoles.h"
#include "reconstruction/eucm_epipolar.h"

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

/*
dynamic algorithm for descriptor comparison

returns a verctor of costs for each possible disparity value
*/
vector<int> compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec, int flawCost);

class EnhancedStereo
{
public:
    
    EnhancedStereo(const EnhancedCamera * cam1, const EnhancedCamera * cam2,
            const StereoParameters & parameters) :
            // initialize members
            params(parameters),
            camera1(cam1->clone()),
            camera2(cam2->clone()),
            epipolarCurves(cam1, cam2, 2000, params.verbosity)
    { 
    }
    
    virtual ~EnhancedStereo()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
    }
    
    
    void setTransformation(const Transf & T12) 
    {
        transf12 = T12; 
        triangulator.setTransformation(T12);
        epipolarCurves.setTransformation(T12);
    }
    
    const Transf & transf() const { return transf12; }
    
//    bool triangulate(double u1, double v1, double u2, double v2, Vector3d & X) const;
    
    double triangulate(double u1, double v1, double u2, double v2, CameraIdx camIdx = CAMERA_1) const;
    
    bool triangulate(const double u1, const double v1, const double u21, const double v21,
            const double u22, const double v22, double & d, double & sigma,
            CameraIdx camIdx = CAMERA_1) const;
    
    const StereoEpipoles & epipoles() const { return epipolarCurves.getEpipoles(); }
    
protected:
    
    EnhancedEpipolar epipolarCurves;
    
    StereoParameters params;
    
    EnhancedCamera *camera1, *camera2;
    
private:
    Triangulator triangulator;
    Transf transf12;  // pose of camera 2 wrt camera 1
};

