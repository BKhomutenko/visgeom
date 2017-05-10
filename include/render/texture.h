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

class Texture
{
public:
    Texture(const Mat8u & textureMat);
        
    virtual ~Texture(); 

    uchar sample(const Vector2d & pt, const Matrix2d & basis) const;

    int cols() const { return _textureMat.cols; }
    int rows() const { return _textureMat.rows; }
    
    const Mat8u _textureMat;
    const Grid2D<uchar> _grid;
    const BiCubicInterpolator<Grid2D<uchar>> _interpolator;
};

