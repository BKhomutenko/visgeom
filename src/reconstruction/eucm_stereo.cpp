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
#include "projection/eucm.h"

#include "reconstruction/scale_parameters.h"


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
