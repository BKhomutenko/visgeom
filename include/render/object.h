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

class IObject
{
public:
    virtual bool intersection(const Vector3d & pos, const Vector3d & dir,
                Vector2d & uv, double & depth) const = 0;
    virtual uchar sample(const Vector2d & pt, const Matrix2d & base) const = 0;
    
    virtual ~IObject() {}
    
};

