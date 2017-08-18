/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your opton) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "calibration/corner_detector.h"

//#include "utils/curve_rasterizer.h"

template<typename T>
void setZero(T begin, T end, T it)
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

bool comp(const pair<double, Vector2i> & a, const pair<double, Vector2i> & b)
{
    return a.first < b.first;
}

//FIXME let it be here for a while
class UnionFind {
  
    int _count;
    vector<int> idVec, sizeVec;
  
public:
    // Create an empty union find data structure with N isolated sets.
    UnionFind(int N) :
        _count(N)
    {
        for (int i = 0; i < _count; i++)
        {
            idVec.push_back(i);
            sizeVec.push_back(1);
        }
    }
    
    ~UnionFind() { }

    // Return the id of component corresponding to object p.
    int findRoot(int p) 
    {
        int root = p;
        while (root != idVec[root]) root = idVec[root];
        while (p != root) 
        { 
            int newp = idVec[p]; 
            idVec[p] = root; 
            p = newp; 
        }
        return root;
    }
    
    // Replace sets containing x and y with their union.
    void merge(int x, int y) 
    {
        int i = findRoot(x); 
        int j = findRoot(y); 
        if (i == j) return;
        // make smaller root point to larger one
        _count--;
        if (sizeVec[i] < sizeVec[j]) 
        { 
            idVec[i] = j;
            sizeVec[j] += sizeVec[i]; 
        }
        else 
        { 
            idVec[j] = i; 
            sizeVec[i] += sizeVec[j]; 
        }
    }
    
    // Are objects x and y in the same set?
    bool connected(int x, int y) 
    { 
        return findRoot(x) == findRoot(y); 
    }
    
    // Return the number of disjoint sets.
    int count() { return _count; }
    int size(int p)
    {
        int i = findRoot(p);
        return sizeVec[i];
    }
};


int main(int argc, char** argv) 
{
    const int Nx = 8;
    const int Ny = 5;
    
    Mat8u frame = imread(argv[1], 0);
    
    CornerDetector detector(Nx, Ny, 3, true, true);
    detector.setImage(frame);
    
    Vector2dVec cornerVec;
    detector.detectPattern(cornerVec);
    
    drawPoints(frame, cornerVec);
    imshow("res", frame);
    waitKey();
    return 0;
}

