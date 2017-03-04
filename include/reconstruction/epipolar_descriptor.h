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
#include "utils/filter.h"

#include "utils/curve_rasterizer.h"

class EpipolarDescriptor
{
public:
    EpipolarDescriptor(int length, int waveThresh, const vector<int> & stepVec) :
            LENGTH(length),
            HALF_LENGTH(length / 2),
            WAVE_THRESH(waveThresh*LENGTH),
            samplingStepVec(stepVec) {}   
                
    // return: the sampling step
    int compute(const Mat8u & img1, const CurveRasterizer<int, Polynomial2> & descRasterRef,
                vector<uint8_t> & descVec)
    {
        descVec.resize(LENGTH);
        bool imageBorder = false;
        bool saturated = false;
        for (int step : samplingStepVec)
        {
            CurveRasterizer<int, Polynomial2> descRaster(descRasterRef);
            descRaster.setStep(-step);
            descRaster.steps(-HALF_LENGTH);
            for (int i = 0; i < LENGTH; i++, descRaster.step())
            {
                if (descRaster.v < 0 or descRaster.v >= img1.rows 
                    or descRaster.u < 0 or descRaster.u >= img1.cols)
                {
                    imageBorder = true;
                    break;
                }
                descVec[i] = img1(descRaster.v, descRaster.u);
//                if (descVec[i] == 255)
//                {
//                    saturated = true;
//                    break;
//                }
            }
            if (imageBorder /*or saturated*/) return -1;
            descResp = totalVariation(descVec.begin(), descVec.end(), int(0));
            descResp = (descResp * /*255*/100) / (int(descVec[HALF_LENGTH]) + /*255*/30);
            if (goodResp()) return step;
        }
        return samplingStepVec.back();
    }
    
    int getResp() { return descResp; }
    
    bool goodResp() { return abs(descResp) > WAVE_THRESH; }
    
private:
    int descResp;
    const int LENGTH;
    const int HALF_LENGTH;
    const int WAVE_THRESH;
    vector<int> samplingStepVec;
};

