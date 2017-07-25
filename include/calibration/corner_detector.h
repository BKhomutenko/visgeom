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
    SubpixelCorner(const Mat32f & gradu, const Mat32f & gradv, const Vector2d pt0,
                 const int steps = 7, const double length = 7);
            
    bool Evaluate(const double* parameters,
                        double* cost,
                        double* gradient) const;

    virtual int NumParameters() const { return 5; }

    const Grid2D<float> _graduGrid, _gradvGrid;
    vector<double> stepVec;
    const double _stepLength;
    const Vector2d _prior;
};

double findMinDistance(const Vector2dVec & cornerVec, const int rows, const int cols);

//TODO make flags instead of bool
class CornerDetector
{
public:
    CornerDetector(int Nx, int Ny, int initRadius = 5, bool improveDetection = true, bool debug = false);
    
    virtual ~CornerDetector() {}
    
    //used in heapsort
    
    
    void setImage(const Mat8u & img);
    
    bool detectPattern(Vector2dVec & ptVec);
    
    //out : _src, _imgrad, _resp, _avgVal
    void computeResponse();
    
    //out : _hypHeap, _detected
    void selectCandidates();
    
    //TODO rewrite
    // /home/bogdan/projects/data/mapping/fluence_calibration/right/000411.png 
    // /home/bogdan/projects/data/mapping/fluence_calibration/right/000402.png
    bool checkCorner(const Vector2i & pt);
    
    bool scaleInvarient(const Vector2i & pt);
    
    //out : _arcVec, _ptVec
    void constructGraph();
    
    vector<int> selectPattern();
    
    vector<int> extractSequence(int idx0, int idx1);
    
    bool verifyDetection(const vector<int> & idxVec);
    
    void improveCorners(Vector2dVec & pointVec) const;   
    
    void initPoin(const Vector2d & pt, double * data) const;
    
    
    
private:

    const bool DEBUG;
    //DETECTION
    Mat32f _resp, _imgrad;
    Mat32f _gradx, _grady;
    Mat8u _src1, _src2;
    double _avgVal;
    
    //for debug
    Mat8u _detected;
    
    //for the graph construction
    Mat16s _idxMap;
    
    vector<vector<int>> _arcVec;
    vector<Vector2i> _ptVec;
    
    //detected strongest maxima sorted by (u + v)
    vector<pair<double, Vector2i>> _hypHeap;
    
    const int _Nx, _Ny;
    const int MAX_CANDIDATE_COUNT;
    
    Mat8u _img;
    
    const bool IMPROVE_DETECTION;
    const int INIT_RADIUS;
    
    //STRUCTURES AND FUNCTIONS
    
    struct TimePoint
    {
        int t;
        int idx;
        int u, v;
        
        TimePoint(int time, int idx, int u, int v): 
            t(time), idx(idx), u(u), v(v) {}
        
        bool friend operator < (const TimePoint & a, const TimePoint & b) 
        {
            return a.t > b.t;
        }
    };
    
    static bool comp(const pair<double, Vector2i> & a, const pair<double, Vector2i> & b)
    {
        return a.first < b.first;
    }
    
    template<typename T>
    static void setZero(T begin, T end, T it)
    {
        auto ref = *it;
        
        T itDir = it;
        do {
            *itDir = 0;
            itDir++;
            if (itDir == end) itDir = begin;
        } while (*itDir * ref > 0);
        
        T itInv = it;
        do {
            *itInv = 0;
            if (itInv == begin) itInv = end;
            itInv--;
        } while (*itInv * ref > 0);
    }
};
