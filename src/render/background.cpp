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

#include "render/background.h"

Background::Background(const ptree & params):
    _texture( imread(params.get<string>("image_name"), 0) )
{
    _u0 = double(_texture.cols()) / 2;
    _v0 = double(_texture.rows()) / 2;
    _fu = _texture.cols() / (2 * M_PI);
    _fv = _texture.rows() / (2 * M_PI);
    _ex << 1, 0, 0;
    _ey << 0, 1, 0;
    _ez << 0, 0, 1;
}

Background::~Background() {}

bool Background::intersection(const Vector3d & pos, const Vector3d & dir,
            Vector2d & uv, double & depth) const
{
    // the mapping from direction to the image is done via yaw-pitch (phi-theta)
    
    // unnormalized sin(phi) and  cos(phi)
    double sphi = _ex.dot(dir); 
    double cphi = _ez.dot(dir); 
    
    // unnormalized sin(theta) and  cos(theta)
    double cth = sqrt(sphi * sphi + cphi * cphi);
    double sth = _ey.dot(dir);
    
    depth = 150; //FIXME
    double phi = atan2(sphi, cphi);
    double th = atan2(sth, cth);
    
    uv[0] = _u0 + _fu * phi;
    uv[1] = _v0 + _fv * th;
    if (uv[0] < 0 or uv[0] > _texture.cols() - 1 
        or uv[1] < 0 or uv[1] > _texture.rows() - 1) return false;
    else return true;
}

uchar Background::sample(const Vector2d & pt, const Matrix2d & base) const
{
    return _texture.sample(pt, base);
}





