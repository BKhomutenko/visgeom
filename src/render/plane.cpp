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

#include "render/plane.h"

Plane::Plane(const Transf & xi, const Mat8u & img):
    _texture(img)
{
    _t = xi.trans();
    Matrix3d R = xi.rotMat();
    _ex = R.col(0);
    _ey = R.col(1);
    _ez = R.col(2);
}

Plane::~Plane() {}

bool Plane::intersection(const Vector3d & pos, const Vector3d & dir,
            Vector2d & uv, double & depth) const
{
    
//        cout << pos.transpose() << "   " << dir.transpose() << endl;
    
    double lambdaNum = _ez.dot(_t - pos);
    double lambdaDen = _ez.dot(dir);
    
    if (lambdaDen < 0.01 ) //TODO condition to rewise, perhaps define the literal elswhere
    {
        return false;
    }
    double lambda = lambdaNum / lambdaDen; 
    depth = lambda * dir.norm();
    Vector3d pt = dir * lambda + pos - _t;
    uv[0] = _u0 + _fu * _ex.dot(pt);
    uv[1] = _v0 + _fv * _ey.dot(pt);
    if (uv[0] < 0 or uv[0] > _texture.cols() - 1 
        or uv[1] < 0 or uv[1] > _texture.rows() - 1) return false;
    else return true;
}

uchar Plane::sample(const Vector2d & pt, const Matrix2d & base) const
{
    return _texture.sample(pt, base);
}




