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

#include "json.h"
#include "ocv.h"
#include "eigen.h"
#include "geometry/geometry.h"

#include "render/object.h"
#include "render/texture.h"

class Plane : public IObject
{
public:

    Plane(const ptree & params);
    
    virtual ~Plane();
    
    bool intersection(const Vector3d & pos, const Vector3d & dir,
                Vector2d & uv, double & depth) const;
    
    uchar sample(const Vector2d & pt, const Matrix2d & base) const;
    
//protected:
    double _u0, _v0, _fu, _fv; // affine transformation to define the origin on the image
    Texture _texture;
    Vector3d _ex, _ey, _ez, _t; //transformation defined by basis vectors and translation
};


