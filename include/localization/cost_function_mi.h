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
Cost functions for localization based on photometric data and mutual information
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"
#include "ceres.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "localization/local_cost_functions.h"
/*
Mutual Information cost function 
//TODO complete the gradient computation explanation
*/
struct MutualInformation : public FirstOrderFunction
{

    MutualInformation(const ICamera * camera, const PhotometricPack & dataPack, const Transf xiBaseCam,
            const Mat32f & img2, double scale, int numBins, double valMax = 1.) :
            _camera(camera->clone()),
            _dataPack(dataPack),
            _xiBaseCam(xiBaseCam),
            _imageGrid(img2.cols, img2.rows, (float*)img2.data),
            _invScale(1. / scale),
            _numBins(numBins),
            _histStep(valMax / (numBins - 1)),
            _increment(1. / dataPack.cloud.size()),
            _hist1(computeHist(dataPack.valVec))
    { }
    
    virtual int NumParameters() const { return 6; }
    
    virtual ~MutualInformation()  
    {
        delete _camera;
        _camera = NULL;
    }
    
    virtual bool Evaluate(double const * parameters, double * cost, double * gradient) const;
    
    void computeShareDerivative(double val, int & idx1, int & idx2, double & der) const;
    
    void computeShares(double val, int & idx1, int & idx2, double & share) const;
    
    vector<double> computeHist(const vector<double> & valVec) const;
    
    //the first vector corresponds to the first image
    //to efficiently compute the histogram for the second image
    vector<double> computeHist2d(const vector<double> & valVec1, const vector<double> & valVec2) const;
    
    vector<double> reduceHist(const vector<double> & hist2d) const;
    
    ICamera * _camera;
    const PhotometricPack & _dataPack;
    const Grid2D<float> _imageGrid;
//    const double _scale;
    const double _invScale;
    
    //histogram params
    const int _numBins;
    const double _histStep;
    
    const Transf _xiBaseCam;
    
    //variable must be initialize for MutualInformation::computeShares
    double _increment;
    
    vector<double> _hist1;
};

struct MutualInformationOdom : public MutualInformation
{

    MutualInformationOdom(  const ICamera * camera, 
                        const PhotometricPack & dataPack, 
                        const Transf xiBaseCam,
                        const Transf xiOdom,
                        const Transf xiPrior,
                        const Mat32f & img2, 
                        double scale, 
                        int numBins, 
                        double valMax = 1.,
                        const double errV = 0.1, 
                        const double errW = 0.01,
                        const double lambdaT = 0.01, 
                        const double lambdaR = 0.01
            );
    
    virtual int NumParameters() const { return 6; }
    
    virtual ~MutualInformationOdom()  
    {
    }
    
    virtual bool Evaluate(double const * parameters, double * cost, double * gradient) const;
    
    Matrix6d _C;
    Matrix6d _J;
    Transf _xiPrior;
};


