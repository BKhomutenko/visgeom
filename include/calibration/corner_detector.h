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
    SubpixelCorner(const Mat32f & gradu, const Mat32f & gradv, const int steps = 7, const double length = 7) :    
            _graduGrid(gradu.cols, gradu.rows, (float*)gradu.data),
            _gradvGrid(gradv.cols, gradv.rows, (float*)gradv.data)
    {
        double stepLength = length / steps;
        stepVec.reserve(2*steps - 2);
        for (int i = 1; i <= steps; i++)
        {
            stepVec.push_back(-i * stepLength);
            stepVec.push_back(i * stepLength);
        }
    }
            
    virtual bool Evaluate(const double* parameters,
                        double* cost,
                        double* gradient) const 
    {
        const double & u = parameters[0];
        const double & v = parameters[1];
        
        //TODO possibly merge
        *cost = 0;
        if (gradient != NULL) fill(gradient, gradient + 4, 0.);
        
        ceres::BiCubicInterpolator<Grid2D> graduInter(_graduGrid);
        ceres::BiCubicInterpolator<Grid2D> gradvInter(_gradvGrid);
        for (int direction = 0; direction < 2; direction++)
        {
            int thIdx = 2 + direction;
            const double s = sin(parameters[thIdx]);
            const double c = cos(parameters[thIdx]);
            double flowDir = direction ? 1 : -1;
            for (int length : stepVec)
            {
                double eta = (length > 0 ? 1 : -1) * flowDir;
                double ui = u + c * length;
                double vi = v + s * length;
                
                double fu, fuu, fuv;
                graduInter.Evaluate(vi , ui,
                            &fu, &fuv, &fuu);
                double fv, fvu, fvv;
                gradvInter.Evaluate(vi , ui,
                            &fv, &fvv, &fvu);
////                cout << fvu << "  " << fuv << endl;            
                *cost += eta*(fu * s - fv * c);
                if (gradient != NULL)
                {
                    gradient[0] += eta * (fuu * s - fvu * c);
                    gradient[1] += eta * (fuv * s - fvv * c);
                    gradient[thIdx] += eta * (-fuu * length * s * s - fvv * length * c * c  
                                    + (fuv + fvu) * length * s * c
                                    + fu * c + fv * s);
                }
            }
        }
        return true;
    }

    virtual int NumParameters() const { return 4; }

    const Grid2D _graduGrid, _gradvGrid;
    vector<double> stepVec;
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
    Mat32f _img;
//    Mat8u _img;
    
    const int INIT_RADIUS;
    
};
