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
Stereo vision definition.
*/

//Standard libraries
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
//Eigen
#include <Eigen/Eigen>

#include "geometry.h"
#include "vision.h"

using namespace std;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;

void computeEssentialMatrix(const vector<Vector3d> & xVec1,
        const vector<Vector3d> & xVec2,
        Matrix3d & E)
{
    assert(xVec1.size() == xVec2.size());
    Matrix<double, 8, 8> AtA = Matrix<double, 8, 8>::Zero();
    Matrix<double, 8, 1> AtB = Matrix<double, 8, 1>::Zero();
    
    for (unsigned int i = 0; i < xVec1.size(); i++)
    {
        const double & x = xVec1[i][0];
        const double & y = xVec1[i][1];
        const double & z = xVec1[i][2];
        const double & u = xVec2[i][0];
        const double & v = xVec2[i][1];
        const double & w = xVec2[i][2];
        Matrix<double, 1, 8> A;
        A << x*u, x*v, x*w, y*u, y*v, y*w, z*u, z*v;

        double b = -z*w;
        AtA += A.transpose() * A;
        AtB += A.transpose() * b;
        if (AtA != AtA)
        {
            cout << i << endl << A << endl;
        }
    }
    
    cout << "E computations" << endl;
    cout << AtA << endl << endl;
    cout << AtB << endl;
    
    Matrix<double, 8, 1> H = AtA.inverse() * AtB;
    E << H(0), H(1), H(2), H(3), H(4), H(5), H(6), H(7), 1;
    cout << E << endl;
}


void StereoSystem::projectPointCloud(const vector<Vector3d> & src,
        vector<Vector2d> & dst1, vector<Vector2d> & dst2) const
{
    dst1.resize(src.size());
    dst2.resize(src.size());
    
    vector<Vector3d> Xc;
    
    TbaseCam1.inverseTransform(src, Xc);
    cam1->projectPointCloud(Xc, dst1);
    
    TbaseCam2.inverseTransform(src, Xc);    
    cam2->projectPointCloud(Xc, dst2);
}

//TODO not finished
bool StereoSystem::triangulate(const Vector3d & v1, const Vector3d & v2,
        const Vector3d & t, Vector3d & X)
{
    //Vector3d v1n = v1 / v1.norm(), v2n = v2 / v2.norm();
    double v1v2 = v1.dot(v2);
    double v1v1 = v1.dot(v1);
    double v2v2 = v2.dot(v2);
    double tv1 = t.dot(v1);
    double tv2 = t.dot(v2);
    double delta = -v1v1 * v2v2 + v1v2 * v1v2;
    if (abs(delta) < 1e-4) // TODO the constant to be revised
    {
        X << -1, -1, -1;
        return false;
    }
    double l1 = (-tv1 * v2v2 + tv2 * v1v2)/delta;
    double l2 = (tv2 * v1v1 - tv1 * v1v2)/delta;
    X = (v1*l1 + t + v2*l2)*0.5;
    return true;
}
     
void StereoSystem::reconstructPointCloud(const vector<Vector2d> & src1,
        const vector<Vector2d> & src2, vector<Vector3d> & dst) const
{
    assert(src1.size() == src2.size());
    dst.resize(src1.size());    
    
    vector<Vector3d> vVec1, vVec2;
    cam1->reconstructPointCloud(src1, vVec1);
    cam2->reconstructPointCloud(src2, vVec2);
    
    TbaseCam1.rotate(vVec1, vVec1);
    TbaseCam2.rotate(vVec2, vVec2);
    Vector3d t = TbaseCam2.trans() - TbaseCam1.trans();
    for (unsigned int i = 0; i < src1.size(); i++)
    {
        Vector3d & X = dst[i];
        triangulate(vVec1[i], vVec2[i], t, X);
        X += TbaseCam1.trans();
    }
}

StereoSystem::~StereoSystem()
{
    if (cam1 != NULL) delete cam1;
    if (cam2 != NULL) delete cam2;
}


