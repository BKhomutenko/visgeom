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

class PlanarGenerator
{
public:
    PlanarGenerator(const ICamera * camera);
    virtual ~PlanarGenerator() { delete _camera; }

    
    virtual uchar brightness(const Vector2d p, const Vector2d u, const Vector2d v) const { return 255; }
    virtual bool isInside(const Vector2d p) const { return false; }
    
    void generate(Mat8u & dst, const Transf xiCam);
    void generateDepth(Mat32f & dst, const Transf xiCam);
    
    void computeDirections();
    
    void computeLengths();
    void computeIntersectionCoordinates();
    void setBackground(const Mat8u & back) { back.copyTo(_back); }
    void setPlaneTransform(const Transf & xiPlane) { _xiPlane = xiPlane; }
    
private:
    ICamera * _camera;
    Transf _xiPlane;
    Vector3dVec _directionVec;
    Mat8u _back;
    
    //temporary variables
    Transf _xiCamPlane;
    Matrix3d _R;
    Vector3d _ex, _ey, _ez, _t;
    double _lambdaNum;
    vector<double> _lambdaVec;
    Vector2dVec _intersectionVec;
    vector<bool> _intersectionFlagVec;
};

class BoardGenerator : public PlanarGenerator
{
public:
    BoardGenerator(const ICamera * camera, int Nx, int Ny, double cellSize);
    virtual ~BoardGenerator() {}

    uchar brightness(const Vector2d p, const Vector2d u, const Vector2d v) const;
    
private:
    int _Nx, _Ny;
    double _cellSize;
};

class ImageGenerator : public PlanarGenerator
{
public:
    ImageGenerator(const ICamera * camera, const Mat8u & img, double scale = 100);
    virtual ~ImageGenerator() {}

    uchar brightness(const Vector2d p, const Vector2d u, const Vector2d v) const;
    
    bool isInside(const Vector2d p) const;
    
private:
    double _scale;
    Vector2d _p0;
    Mat8u _img;
    Mat32s _J;
};
