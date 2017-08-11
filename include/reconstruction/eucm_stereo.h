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
#include "json.h"
#include "geometry/geometry.h"
#include "projection/eucm.h"

#include "reconstruction/scale_parameters.h"
#include "reconstruction/triangulator.h"
#include "reconstruction/epipoles.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/eucm_epipolar.h"

struct StereoParameters : public ScaleParameters
{
//    StereoParameters() {}
    StereoParameters(const ScaleParameters & scaleParams) : ScaleParameters(scaleParams) {}
    
    StereoParameters(const ptree & params);
    
    int dispMax = 48;
  
    int maxError = 25;
//    int maxBias = 10;
    
    int verbosity = 0;

    //squared value, corresonds to 50 px
    int epipoleMargin = 2500;

    //multi-hypoteses
    //TODO remove?
    int hypMax = 1;
    int maxHypDiff = 10;
    
    //for descriptor matching, depends on the image noise
    int flawCost = 7;
    
    //descriptor computation parameters
    int descLength = 5;
    vector<int> scaleVec = {1, 2, 3, 5};
    int descRespThresh = 3;
    
    int numEpipolarPlanes = 2000;
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
            const StereoParameters & params);
    
    virtual ~EnhancedStereo();
    
    void setTransformation(const Transf & T12);    
    const Transf & transf() const { return _transf12; }
    const Matrix3d & R12() const { return _R12; }
    const Matrix3d & R21() const { return _R21; }
    const Vector3d & t12() const { return _t12; }
    
//    bool triangulate(double u1, double v1, double u2, double v2, Vector3d & X) const;
    
    double triangulate(double u1, double v1, double u2, double v2, CameraIdx camIdx = CAMERA_1) const;
    
    //TODO remove this function, obsolete 
    bool triangulate(const double u1, const double v1, const double u21, const double v21,
            const double u22, const double v22, double & d, double & sigma,
            CameraIdx camIdx = CAMERA_1) const;
    
    bool triangulate(const Vector2d pt11, const Vector2d pt12, const Vector2d pt21,
        const Vector2d pt22, double & d, double & sigma) const;
        
    const StereoEpipoles & epipoles() const { return _epipolarCurves.getEpipoles(); }
    
 
protected:
    StereoParameters _params;
    EnhancedCamera *_camera1, *_camera2;
    EnhancedEpipolar _epipolarCurves;
    EpipolarDescriptor _epipolarDescriptor;
    
    const int HALF_LENGTH;
    const int MARGIN;
    
private:
    Triangulator _triangulator;
    Transf _transf12;  // pose of camera 2 wrt camera 1
    Matrix3d _R12, _R21;
    Vector3d _t12;    
};

