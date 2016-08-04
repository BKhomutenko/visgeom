#pragma once

#include "eigen.h"


// Multi-hypothesis point-cloud data structure
struct MHPack
{
    std::vector<int> idxVec; // Index of the point in the depthmap vector
    Vector2dVec imagePointVec; // Index of the pixel on the original image
    std::vector<int> hypIdxVec; // Index for the hypothesis from the original depthmap
    std::vector<double> costVec; // Cost of each hypothesis
    Vector3dVec cloud; // Cloud points - exact data specified by datatype flag
    std::vector<double> valVec; // Scalar values - exact data specified by datatype flag

    enum Data
    {
    	RECONSTRUCTION_WITH_IMAGE_VALUES,
    	RECONSTRUCTION_WITH_SIGMA,
    	MINMAX_DISTANCE_VEC_WITH_EMPTY
    };

    Data datatype; // datatype of the values in cloud and valVec
};