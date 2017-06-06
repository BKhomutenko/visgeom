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
    
#pragma once

#include "eigen.h"
#include "ceres.h"
#include "ocv.h"
//#include "io.h"

class SubpixelCorner : public ceres::FirstOrderFunction 
{
public:
    SubpixelCorner(const Mat32f & gradu, const Mat32f & gradv, const int steps = 7, const double length = 7);
            
    bool Evaluate(const double* parameters,
                        double* cost,
                        double* gradient) const;

    virtual int NumParameters() const { return 5; }

    const Grid2D<float> _graduGrid, _gradvGrid;
    vector<double> stepVec;
    const double _stepLength;
};

double findMinDistance(const Vector2dVec & cornerVec, const int rows, const int cols);

class CornerDetector
{
public:
    CornerDetector(const Mat8u & img, const int initRadius = 5);
    
    virtual ~CornerDetector() {}
    
    void improveCorners(Vector2dVec & pointVec) const;   
    
    void initPoin(const Vector2d & pt, double * data) const;
    
private:
    Mat32f _gradx, _grady;
//    Mat32f _img;
    Mat8u _img;
    
    const int INIT_RADIUS;
    
};
