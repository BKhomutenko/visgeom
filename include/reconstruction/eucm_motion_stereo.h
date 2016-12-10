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
Depth-from-motion class for semidense depth estimation
*/

#pragma once


#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "utils/filter.h"
#include "geometry/geometry.h"
#include "projection/eucm.h"

#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_stereo.h"

//TODO add errorMax threshold
struct MotionStereoParameters : StereoParameters
{
    int descLength = 5;
    int gradientThresh = 3;
};


class MotionStereo : private EnhancedStereo
{
public:
    MotionStereo(const EnhancedCamera * cam1, 
        const EnhancedCamera * cam2, MotionStereoParameters parameters) :
        EnhancedStereo(cam1, cam2, parameters),
        params(parameters)
    {
    }

    virtual ~MotionStereo()
    {    
    }    
    
    //TODO figure out how to treat the mask efficiently
    void setBaseImage(const Mat8u & image)
    {
        image.copyTo(img1);
        computeMask();
    }
       
    /*
    -Select salient points and points with defined depth
    -reproject all the points onto the next image
    -for those which have small uncertainties just keep that value
    -for those with wide uncertainty or without a value recompute it

    hypothesis quality must be evaluated using normal Kalman filtering
    That is, at every step the uncertainty grows. Once it reaches a certain value,
    the hypothesis is removed

    two sets of points must be treated separately:
    -points with bad descriptors but with depth estimation. 
        Project forward, increment uncertainty
    -points with good descriptors : 
        -with depth estimation. Project forward using serach 
        (if the search distance is greater than 1)
        -withoud depth estimation. Generate new hypotheses using complete epipolar search
    */
    void reprojectDepth(Transf T12, const Mat8u & img2, DepthMap & depth);
    
    //FIXME for backward compatibility, to be fixed or removed
    void computeDepth(Transf T12, const Mat8u & img2, DepthMap & depth);
    
    void validateDepth(Transf T12, const Mat8u & img2, DepthMap & depth);
       
private:
    
    // based on the image gradient
    void computeMask()
    {
        Mat16s gradx, grady;
        Sobel(img1, gradx, CV_16S, 1, 0, 1);
        Sobel(img1, grady, CV_16S, 0, 1, 1);
        Mat16s gradAbs = abs(gradx) + abs(grady);
        GaussianBlur(gradAbs, gradAbs, Size(7, 7), 0, 0);
        Mat8u gradAbs8u;
        gradAbs.convertTo(gradAbs8u, CV_8U);
        threshold(gradAbs8u, maskMat, params.gradientThresh, 128, CV_THRESH_BINARY);
    }
   
    Mat8u img1;    
    Mat8u maskMat; //TODO compute mask
    const MotionStereoParameters params;
};

