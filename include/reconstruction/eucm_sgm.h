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
#include "ocv.h"
#include "eigen.h"
#include "json.h"

#include "geometry/geometry.h"
#include "projection/eucm.h"

#include "reconstruction/scale_parameters.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_stereo.h"

struct SgmParameters : public StereoParameters
{
    SgmParameters(const ptree & params):
            StereoParameters(params)
    {
        for (auto & item : params.get_child("sgm_stereo_parameters"))
        {
            const string & pname = item.first;
            if (pname == "step_cost")               lambdaStep = item.second.get_value<int>();
            else if (pname == "jump_cost")          lambdaJump = item.second.get_value<int>();
            else if (pname == "image_based_cost")       imageBasedCost = item.second.get_value<bool>();
            else if (pname == "salient_points_only")    salientPoints = item.second.get_value<bool>();
            else if (pname == "use_uv_cache")           useUVCache = item.second.get_value<bool>();
        }
    }
    
    // cost parameters
    int lambdaStep = 5;
    int lambdaJump = 32;
    
    bool imageBasedCost = true;
    
    //TODO get a meaningful threshold/ method to discard
    //temporarily descriptor step must be 1
    //all non-salient points are discarded
    bool salientPoints = true;
    
    //precompute all the epipolar curves
    bool useUVCache = true;
};

//TODO revamp, take MotionStereo as a model
class EnhancedSgm : private EnhancedStereo
{
public:
    
    EnhancedSgm(Transf T12, const EnhancedCamera * cam1,
            const EnhancedCamera * cam2, const SgmParameters & params) :
            EnhancedStereo(cam1, cam2, params),
            _params(params)
    { 
        setTransformation(T12);
        assert(params.dispMax % 2 == 0);
        createBuffer();
        computeReconstructed();
        computeRotated();
        computePinf();
        if (params.useUVCache) computeUVCache();
    }
    
    virtual ~EnhancedSgm()
    {
    }
    
    // precompute coordinates for different disparities to speedup the computation
    void computeUVCache();
    
    // An interface function
    void computeStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depthMap);
    
    //// EPIPOLAR GEOMETRY
    
    // computes reconstVec -- reconstruction of every pixel of the first image
    void computeReconstructed();
    
    // computes reconstRotVec -- reconstVec rotated into the second frame
    void computeRotated();
       
    // computes pinfVec -- projections of all the reconstructed points from the first image
    // onto the second image as if they were at infinity
    void computePinf();
    
    // calculate the coefficients of the polynomials for all the 
    void computeEpipolarIndices();
    
    //// DYNAMIC PROGRAMMING
    void createBuffer();
       
    // fill up the error buffer using 2*S-1 pixs along epipolar lines as local desctiprtors
    void computeCurveCost(const Mat8u & img1, const Mat8u & img2);
    
    void computeDynamicProgramming();
    
    void computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost);
    
    void reconstructDisparityMH();
    void reconstructDisparity();  // using the result of the dynamic programming
    
    //TODO rewrite
    void reconstructDepth(DepthMap & depth) const;
    //// MISCELLANEOUS
    
    // index of an object in a linear array corresponding to pixel [row, col] 
    int getLinearIndex(int x, int y) const { return _params.xMax*y + x; }
      
    CurveRasterizer<int, Polynomial2> getCurveRasteriser(CameraIdx camIdx, int idx,
                                                         uint32_t * flags = NULL) const;
    
    void skipPixel(int x, int y);
    
    // reconstruction
    void fillGaps(uint8_t * const data, const int step);
    
    Mat32s & disparity() { return _smallDisparity; }
    
private:
    
    std::vector<bool> _maskVec;
    
    Vector2dVec _pointVec1;  // the depth points on the image 1
    Vector3dVec _reconstVec;  // reconstruction of every pixel by cam1
    Vector3dVec _reconstRotVec;  // reconstVec rotated into the second frame
    Vector2dVec _pinfVec;  // projection of reconstRotVec by cam2
    
    // discretized version
    Vector2iVec _pointPxVec1;
    Vector2iVec _pinfPxVec;
    
    // to be able to change the jump cost
    int _jumpCost;
    
    const int DISPARITY_MARGIN = 20;
    Mat32s _uCache, _vCache;
    Mat8u _errorBuffer;
    Mat8u _costBuffer; //TODO maybe merge with salientBuffer
    Mat8u _salientBuffer; 
    Mat8u _stepBuffer;
    Mat8u _skipBuffer;
    Mat32s _tableauLeft, _tableauRight; //FIXME check the type through the code
    Mat32s _tableauTop, _tableauBottom;
    Mat32s _smallDisparity;
    Mat32s _finalErrorMat;
    
    
    const SgmParameters _params;
};

