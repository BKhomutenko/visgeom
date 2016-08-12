#pragma once

#include "eigen.h"


// Multi-hypothesis point-cloud data structure
struct MHPack
{
    std::vector<int> idxVec; // Index of the point in the depthmap vector
    Vector2dVec imagePointVec; // Index of the pixel on the original image
    std::vector<int> hypIdxVec; // Index for the hypothesis from the original depthmap
    std::vector<double> costVec; // Cost of each hypothesis
    std::vector<double> sigmaVec; // Position uncertainty of each hypothesis
    Vector3dVec cloud; // Cloud points - exact data specified by datatype flag
    std::vector<double> valVec; // Image values at each position

    enum Data // TODO - Fix or remove this, as this is no longer necessary except for specifying MIXMAX
    {
    	NO_DATA,
        RECONSTRUCTION_WITH_IMAGE_VALUES,
    	RECONSTRUCTION_WITH_SIGMA,
    	MINMAX_DISTANCE_VEC_WITH_SIGN  // In this case, `cloud` will have twice the number of points, in nearest-farthest order
    };

    Data datatype; // datatype of the values in cloud

    // Default constructor
    MHPack() :
        datatype(NO_DATA) {}

    // Copy constructor
    MHPack(const MHPack & other) :
        idxVec(other.idxVec),
        imagePointVec(other.imagePointVec),
        hypIdxVec(other.hypIdxVec),
        costVec(other.costVec),
        sigmaVec(other.sigmaVec),
        cloud(other.cloud),
        valVec(other.valVec),
        datatype(other.datatype) {}
};