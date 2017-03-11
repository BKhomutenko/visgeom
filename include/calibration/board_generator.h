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

#include "std.h"
#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "projection/generic_camera.h"

class BoardGenerator
{
public:
    BoardGenerator(const ICamera * camera, int Nx, int Ny, double cellSize, double factor = 10);
    virtual ~BoardGenerator() { delete _camera; }

    void generate(Mat8u & dst, const Transf xi) const;
    
    uchar computeBrightness(const double x, const double y, const double smooth) const;
    
private:
    ICamera * _camera;
    int _Nx, _Ny;
    double _cellSize;
};
