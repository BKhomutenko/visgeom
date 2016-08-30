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
#include "ocv.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "camera/eucm.h"

#include "utils/scale_parameters.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/epipoles.h"


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
    double delta = -pp * qq + pq * pq;
    if (abs(delta) < 1e-10) // TODO the constant to be revised
    {
        if (params.verbosity > 2) 
        {
            cout << "    not triangulated " << abs(delta) << " " << (abs(delta) < 1e-10) << endl;
        }
        X = Vector3d(0, 0, 0);
        return false;
    }
    double l1 = (-tp * qq + tq * pq)/delta;
    double l2 = (tq * pp - tp * pq)/delta;
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
        return -1;
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
    double delta = -pp * qq + pq * pq;
    if (abs(delta) < 1e-10) // TODO the constant to be revised
    {
        if (params.verbosity > 2) 
        {
            cout << "    not triangulated " << abs(delta) << " " << (abs(delta) < 1e-10) << endl;
        }
        return -1;
    }
    if (camIdx == CAMERA_1) return (-tp * qq + tq * pq)/delta;
    else return (tq * pp - tp * pq)/delta;
}


vector<int> EnhancedStereo::compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec) const
{
    const int HALF_LENGTH = desc.size() / 2;
    vector<int> rowA(sampleVec.size()), rowB(sampleVec.size());
    
    //match the first half
    for (int i = 0; i < sampleVec.size(); i++)
    {
        rowA[i] = (abs(sampleVec[i] - desc[0]));
    }
    for (int i = 1; i <= HALF_LENGTH; i++)
    {
        rowB[0] = (rowA[0] + params.flawCost + abs(sampleVec[0] - desc[i]));
        int cost = min(rowA[i] + params.flawCost, rowA[0]);
        rowB[1] = (cost + abs(sampleVec[i] - desc[i]));
        for (int j = 2; j < sampleVec.size(); j++)
        {
            cost = min(min(rowA[j] + params.flawCost, rowA[j - 1]), rowA[j - 2] + params.flawCost);
            rowB[j] = (cost+ abs(sampleVec[j] - desc[i]));
        }
        swap(rowA, rowB);
    }
    vector<int> rowC(sampleVec.size()); //center cost
    swap(rowA, rowC);
    
    //match the second half (from the last pixel to first)
    for (int i = 0; i < sampleVec.size(); i++)
    {
        rowA[i] = (abs(sampleVec[i] - desc.back()));
    }
    for (int i = desc.size() - 1; i > HALF_LENGTH + 1; i--)
    {
        for (int j = 0; j < sampleVec.size() - 2; j++)
        {
            int cost = min(min(rowA[j] + params.flawCost, rowA[j + 1]), rowA[j + 2] + params.flawCost);
            rowB[j] = (cost + abs(sampleVec[j] - desc[i]));
        }
        int j = sampleVec.size() - 2;
        int cost = min(rowA[j] + params.flawCost, rowA[j + 1]);
        rowB[j] = (cost + abs(sampleVec[j] - desc[i]));
        rowB.back() = (rowA.back() + params.flawCost + abs(sampleVec.back() - desc[i]));
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
    
