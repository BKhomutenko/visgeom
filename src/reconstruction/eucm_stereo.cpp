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

#include "reconstruction/eucm_stereo.h"

#include "std.h"
#include "io.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "camera/eucm.h"

#include "utils/scale_parameters.h"


bool EnhancedStereo::triangulate(double u1, double v1, double u2, double v2, Vector3d & X) const
{
    if (params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    Vector3d p, q;
    if (not camera1->reconstructPoint(Vector2d(u1, v1), p) or 
        not camera2->reconstructPoint(Vector2d(u2, v2), q) )
    {
        if (params.verbosity > 2) 
        {
            cout << "    not reconstructed " << Vector2d(u1, v1).transpose(); 
            cout << " # " << Vector2d(u2, v2).transpose() << endl;
        }
        X = Vector3d(0, 0, 0);
        return false;
    }
    Vector3d t = Transform12.trans();
    q = Transform12.rotMat() * q;
    if (params.verbosity > 3) 
    {
        cout << "    pt1: " << u1 << " " << v1 << endl;
        cout << "    p: " << p.transpose() << endl;
        cout << "    pt2: " << u2 << " " << v2 << endl;
        cout << "    q: " << q.transpose() << endl;
    }
    double pq = p.dot(q);
    double pp = p.dot(p);
    double qq = q.dot(q);
    double tp = t.dot(p);
    double tq = t.dot(q);
    double ppqq = pp * qq;
    double delta = -pp * qq + pq * pq;
    if (abs(delta) < ppqq * 1e-3) // TODO the constant to be revised
    {
        X = p * MAX_DEPTH;
        return true;
    }
    double deltainv = 1/delta;
    double l1 = (-tp * qq + tq * pq) * deltainv;
    double l2 = (tq * pp - tp * pq) * deltainv;
    X = (p*l1 + t + q*l2)*0.5;
    return true;
}
  
//TODO replace distance by length coefficient
// returns the distance between corresponding camera and the point
double EnhancedStereo::triangulate(double u1, double v1, double u2, double v2,
        CameraIdx camIdx) const
{
    if (params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    Vector3d p, q;
    if (not camera1->reconstructPoint(Vector2d(u1, v1), p) or 
        not camera2->reconstructPoint(Vector2d(u2, v2), q) )
    {
        if (params.verbosity > 2) 
        {
            cout << "    not reconstructed " << Vector2d(u1, v1).transpose(); 
            cout << " # " << Vector2d(u2, v2).transpose() << endl;
        }
        return OUT_OF_RANGE;
    }
    Vector3d t = Transform12.trans();
    q = Transform12.rotMat() * q;
    if (params.verbosity > 4) 
    {
        cout << "    pt1: " << u1 << " " << v1 << endl;
        cout << "    p: " << p.transpose() << endl;
        cout << "    pt2: " << u2 << " " << v2 << endl;
        cout << "    q: " << q.transpose() << endl;
    }
    double pq = p.dot(q);
    double pp = p.dot(p);
    double qq = q.dot(q);
    double tp = t.dot(p);
    double tq = t.dot(q);
    double ppqq = pp * qq;
    double delta = -pp * qq + pq * pq;
    if (abs(delta) < ppqq * 1e-6) // max 100m for 10cm base
    {
        return MAX_DEPTH;
    }
    if (camIdx == CAMERA_1) return sqrt(pp) * (-tp * qq + tq * pq)/delta;
    else return sqrt(qq) * (tq * pp - tp * pq)/delta;
}


int computeError(int v, int thMin, int thMax)
{
    return max(0, max(thMin - v, v - thMax));
}


vector<int> EnhancedStereo::compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec) const
{
    vector<int> thMinVec(sampleVec.size()), thMaxVec(sampleVec.size());
    
    for (int i = 1; i < sampleVec.size() - 1; i++)
    {
        const int & d = desc[i];
        int d1 = (desc[i] + desc[i - 1]) / 2;
        int d2 = (desc[i] + desc[i + 1]) / 2;
        thMinVec[i] = min(d, min(d1, d2));
        thMaxVec[i] = max(d, min(d1, d2));
    }
    
    if (desc[0] > desc[1])
    {
        thMinVec[0] = (desc[0] + desc[1]) / 2;
        thMaxVec[0] = desc[0]; 
    }
    else
    {
        thMaxVec[0] = (desc[0] + desc[1]) / 2;
        thMinVec[0] = desc[0];
    }
    
    const int & d = desc[sampleVec.size() - 2];
    if (desc.back() > d)
    {
        thMinVec.back() = (desc.back() + d) / 2;
        thMaxVec.back() = desc.back(); 
    }
    else
    {
        thMaxVec.back() = (desc.back() + d) / 2;
        thMinVec.back() = desc.back();
    }
    
    
    const int HALF_LENGTH = desc.size() / 2;
    vector<int> rowA(sampleVec.size()), rowB(sampleVec.size());
    
    //match the first half
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(int(sampleVec[i]) - int(desc[0]));
        rowA[i] = computeError(sampleVec[i], thMinVec[0], thMaxVec[0]);
    }
    for (int i = 1; i <= HALF_LENGTH; i++)
    {
//        rowB[0] = rowA[0] + params.flawCost + abs(int(sampleVec[0]) - int(desc[i]));
        rowB[0] = rowA[0] + params.flawCost + computeError(sampleVec[0], thMinVec[i], thMaxVec[i]);
        int cost = min(rowA[i] + params.flawCost, rowA[0]);
//        rowB[1] = cost + abs(int(sampleVec[1]) - int(desc[i]));
        rowB[1] = cost + computeError(sampleVec[1], thMinVec[i], thMaxVec[i]);
        for (int j = 2; j < sampleVec.size(); j++)
        {
            cost = min(min(rowA[j] + params.flawCost, rowA[j - 1]), rowA[j - 2] + params.flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        swap(rowA, rowB);
    }
    vector<int> rowC(sampleVec.size()); //center cost
    swap(rowA, rowC);
    
    //match the second half (from the last pixel to first)
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(sampleVec[i] - desc.back());
        rowA[i] = computeError(sampleVec[i], thMinVec.back(), thMaxVec.back());
    }
    for (int i = desc.size() - 2; i > HALF_LENGTH; i--)
    {
        for (int j = 0; j < sampleVec.size() - 2; j++)
        {
            int cost = min(min(rowA[j] + params.flawCost, rowA[j + 1]), rowA[j + 2] + params.flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        int j = sampleVec.size() - 2;
        int cost = min(rowA[j] + params.flawCost, rowA[j + 1]);
//        rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        
//        rowB.back() = rowA.back() + params.flawCost + abs(int(sampleVec.back()) - int(desc[i]));
        rowB.back() = rowA.back() + params.flawCost + computeError(sampleVec.back(), thMinVec[i], thMaxVec[i]);
        swap(rowA, rowB);
    }
    
    //accumulate the cost
    for (int i = 0; i < sampleVec.size() - 2; i++)
    {
        rowC[i] += min(min(rowA[i] + params.flawCost, rowA[i + 1]), rowA[i + 2] + params.flawCost);
    }
    int i = rowC.size() - 2;
    rowC[i] += min(rowA[i] + params.flawCost, rowA[i + 1]);
    rowC.back() += rowA.back() + params.flawCost;
    return rowC;
} 

