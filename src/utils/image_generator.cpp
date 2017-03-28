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

#include "utils/image_generator.h"

#include "std.h"
#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "projection/eucm.h"

PlanarGenerator::PlanarGenerator(const ICamera * camera) : 
    _camera(camera->clone())
    {
        computeDirections();
    }


void PlanarGenerator::computeDirections()
{
    const int W = _camera->width;
    const int H = _camera->height;
    _directionVec.resize(_camera->width * _camera->height);
    for (int r = 0; r < H; r++)
    {
        const int shift = r * W;
        for (int c = 0; c < W; c++)
        {
            Vector2d pt(c, r);
            if (not _camera->reconstructPoint(pt, _directionVec[shift + c])) 
            {
                _directionVec[shift + c] = Vector3d(0, 0, 0);
            }
        }
    }
}

void PlanarGenerator::computeLengths()
{
    _lambdaVec.resize(_camera->width * _camera->height);
    for (int idx = 0; idx < _directionVec.size(); idx++)
    {
        const double lambdaDenom = _directionVec[idx].dot(_ez);
        if (abs(lambdaDenom) > 0.01 )
        {
            _lambdaVec[idx] = _lambdaNum / lambdaDenom;
        }
        else
        {
            _lambdaVec[idx] = 0;
        }
    }
}

void PlanarGenerator::computeIntersectionCoordinates()
{
    _intersectionVec.resize(_camera->width * _camera->height);
    _intersectionFlagVec.resize(_camera->width * _camera->height);
    fill(_intersectionFlagVec.begin(), _intersectionFlagVec.end(), true);
    for (int idx = 0; idx < _directionVec.size(); idx++)
    {
        if (_lambdaVec[idx] <= 0)
        {
            _intersectionFlagVec[idx] = false;
            continue;
        }
        Vector3d v = _directionVec[idx] * _lambdaVec[idx] - _t;
        _intersectionVec[idx] = Vector2d(v.dot(_ex), v.dot(_ey));
    }
}

void PlanarGenerator::generateDepth(Mat32f & dst, const Transf xiCam)
{
    const int W = _camera->width;
    const int H = _camera->height;
    dst.create(H, W);
    _xiCamPlane = xiCam.inverseCompose(_xiPlane);
    _R = _xiCamPlane.rotMat();
    _ex = _R.col(0);
    _ey = _R.col(1);
    _ez = _R.col(2);
    _t = _xiCamPlane.trans();
    _lambdaNum = _t.dot(_ez);
//    Matrix3d Rcam = xiCam.rotMat();
    computeLengths();
    computeIntersectionCoordinates();
    
    for (int r = 0; r < H; r++)
    {
        const int shift = r * W;
        for (int c = 0; c < W; c++)
        {
            const int idx = shift + c;
            if (_intersectionFlagVec[idx] and isInside(_intersectionVec[idx])) 
            {
                dst(r, c) = _directionVec[idx].norm() * _lambdaVec[idx];
            }
            else
            {
                dst(r, c) = 0;
            }
        }
    }
    
}

void PlanarGenerator::generate(Mat8u & dst, const Transf xiCam)
{
    const int W = _camera->width;
    const int H = _camera->height;
    dst.create(H, W);
    _xiCamPlane = xiCam.inverseCompose(_xiPlane);
    _R = _xiCamPlane.rotMat();
    _ex = _R.col(0);
    _ey = _R.col(1);
    _ez = _R.col(2);
    _t = _xiCamPlane.trans();
    _lambdaNum = _t.dot(_ez);
//    Matrix3d Rcam = xiCam.rotMat();
    computeLengths();
    computeIntersectionCoordinates();
    
//    array<double, 6> backCamParams {0.5, 1, _back.cols / 5, _back.cols / 5, _back.cols / 2, _back.rows / 2};
//    EnhancedCamera backCam(backCamParams.data());
    
    for (int r = 0; r < H; r++)
    {
        cout << r << endl;
        const int shift = r * W;
        for (int c = 0; c < W; c++)
        {
            const int idx = shift + c;
            bool computeBackground = false;
            if (not _intersectionFlagVec[idx]) 
            {
                if (_back.empty()) dst(r, c) = 255;
                else computeBackground = true;
            }
            else
            {
                const Vector2d & p = _intersectionVec[idx];
                Vector2d u, v;
                if (r == 0 or not _intersectionFlagVec[idx - W])
                {
                    v = _intersectionVec[idx + W] - p;
                }
                else if (r == H - 1 or not _intersectionFlagVec[idx + W])
                {
                    v = p - _intersectionVec[idx - W];
                }
                else
                {
                    v = 0.5 * (_intersectionVec[idx + W] - _intersectionVec[idx - W]);
                }
                if (c == 0 or not _intersectionFlagVec[idx - 1])
                {
                    u = _intersectionVec[idx + 1] - p;
                }
                else if (c == W - 1 or not _intersectionFlagVec[idx + 1])
                {
                    u = p - _intersectionVec[idx - 1];
                }
                else
                {
                    u = 0.5 * (_intersectionVec[idx + 1] - _intersectionVec[idx - 1]);
                }
                
                uchar res = brightness(p, u, v);
                if (res == 255)
                {
                    if (_back.empty()) dst(r, c) = 255;
                    else computeBackground = true;
                }
                else dst(r, c) = res;
            }
            
            if (computeBackground)
            {
                throw;                
//                backCam.projectPoint(Rcam * _directionVec[idx]
            }
        }
    }
    //reconstruct
    //find the intersection
    //find the color
}

BoardGenerator::BoardGenerator(const ICamera * camera, int Nx, int Ny,
        double cellSize):
            PlanarGenerator(camera),
            _Nx(Nx),
            _Ny(Ny),
            _cellSize(cellSize) {}

uchar BoardGenerator::brightness(const Vector2d p, const Vector2d u, const Vector2d v) const
{
    double smooth = min(u.norm(), v.norm());
    if (p[0] <= -_cellSize or p[0] >= _Nx * _cellSize) return 255;
    else if (p[1] <= -_cellSize or p[1] >= _Ny * _cellSize) return 255;
    
    const double xrel = p[0] / _cellSize;
    const double yrel = p[1] / _cellSize;
    int cx = ceil(xrel);
    int cy = ceil(yrel);
    
    const double dx = abs(xrel - round(xrel));
    const double dy = abs(yrel - round(yrel));
    const double delta = min(dx, dy);
    if (cx % 2 == cy % 2)
    {
        if (cx != 0 and cy != 0 and cx != _Nx and cy != _Ny)
        {
            const double res = 180*(1 + delta * _cellSize / smooth);
            return uchar(min(res, 255.));
        }
        else return 255;
    }
    else
    {
        const double res = 180*(1 - delta * _cellSize / smooth);
        return uchar(max(res, 0.));
    }
//    double res = 1 - sin(x/_cellSize * M_PI) * sin(y/_cellSize * M_PI)  * _factor;
//    res = max( min(1., (res)), 0. );
//    return uchar(res * 255); 
}


ImageGenerator::ImageGenerator(const ICamera * camera,
    const Mat8u & img,
    double scale):
            PlanarGenerator(camera),
            _scale(scale)
{
    _p0[0] = 0.5 * img.cols / scale;
    _p0[1] = 0.5 * img.rows / scale;
    cv::integral(img, _J);
    img.copyTo(_img);
}

void clip(const vector<Vector2d> & simplex, list<Vector2d> & poly2);

bool ImageGenerator::isInside(const Vector2d p) const
{
    Vector2i p1 = round((p + _p0) * _scale);
    return  (p1[0] >= 0 and p1[0] < _img.cols and p1[1] >= 0 and p1[1] < _img.rows);
}

uchar ImageGenerator::brightness(const Vector2d p, const Vector2d u, const Vector2d v) const
{
    
    //check the image limits
    
    if (not isInside(p)) return 255;
    const double invScale = 1. / _scale;
    //Nneigh
    if (max( max( abs(u[0]), abs(u[1]) ), max( abs(v[0]), abs(v[1]) ) ) < 0.2*invScale )
    {    
        Vector2i p1 = round((p + _p0) * _scale);
        return _img(p1[1], p1[0]);
    }
    
    //box filters
    /*
    array<double, 4> xArr, yArr;
    xArr[0] = p[0] - u[0] * 0.5;
    xArr[1] = p[0] + u[0] * 0.5;
    xArr[2] = p[0] - v[0] * 0.5;
    xArr[3] = p[0] + v[0] * 0.5;
    
    yArr[0] = p[1] - u[1] * 0.5;
    yArr[1] = p[1] + u[1] * 0.5;
    yArr[2] = p[1] - v[1] * 0.5;
    yArr[3] = p[1] + v[1] * 0.5;
    
    double xMin = *min_element(xArr.begin(), xArr.end());
    double yMin = *min_element(yArr.begin(), yArr.end());
    
    double xMax = *max_element(xArr.begin(), xArr.end());
    double yMax = *max_element(yArr.begin(), yArr.end());
    
    int u1 = round((_p0[0] + xMin) * _scale);
    int u2 = round((_p0[0] + xMax) * _scale) + 1;
    int v1 = round((_p0[1] + yMin) * _scale);
    int v2 = round((_p0[1] + yMax) * _scale) + 1;
    
    u1 = min(max(0, u1), _J.cols - 1);
    u2 = min(max(0, u2), _J.cols - 1);
    v1 = min(max(0, v1), _J.rows - 1);
    v2 = min(max(0, v2), _J.rows - 1);
    if (u1 == u2 or v1 == v2) return 255;
    double res = _J(v1, u1) + _J(v2, u2) - _J(v1, u2) - _J(v2, u1);
    return uchar( res / ((u2 - u1) * (v2 - v1)) );
    */
    
    //sinc filter
    array<double, 4> xArr, yArr;
    xArr[0] = p[0] - u[0] * 3;
    xArr[1] = p[0] + u[0] * 3;
    xArr[2] = p[0] - v[0] * 3;
    xArr[3] = p[0] + v[0] * 3;
    
    yArr[0] = p[1] - u[1] * 3;
    yArr[1] = p[1] + u[1] * 3;
    yArr[2] = p[1] - v[1] * 3;
    yArr[3] = p[1] + v[1] * 3;
    
    double xMin = *min_element(xArr.begin(), xArr.end());
    double yMin = *min_element(yArr.begin(), yArr.end());
    
    double xMax = *max_element(xArr.begin(), xArr.end());
    double yMax = *max_element(yArr.begin(), yArr.end());
    
    int u1 = round((_p0[0] + xMin) * _scale);
    int u2 = round((_p0[0] + xMax) * _scale) + 1;
    int v1 = round((_p0[1] + yMin) * _scale);
    int v2 = round((_p0[1] + yMax) * _scale) + 1;
    
    u1 = min(max(0, u1), _img.cols - 1);
    u2 = min(max(0, u2), _img.cols - 1);
    v1 = min(max(0, v1), _img.rows - 1);
    v2 = min(max(0, v2), _img.rows - 1);
    if (u1 == u2 or v1 == v2) return 255;
    
    
    double res = 0;
    double normFactor = 0;
    Matrix2d A;
    A <<    u[0],   v[0],
            u[1],   v[1];
    Matrix2d Ainv = A.inverse();     
    for (int i = v1; i < v2; i++)
    {
        for (int j = u1; j < u2; j++)
        {
            double k = 0;
            for (double di = -0.3; di < 0.51; di+= 0.3)
            {
                for (double dj = -0.3; dj < 0.51; dj+= 0.3)
                {
                    Vector2d delta = Vector2d(j + dj, i + di) * invScale - _p0 - p;
                    Vector2d t = Ainv * delta;
                    k += sinc(M_PI * t[0]) * sinc(M_PI * t[1]);
                }
            }
            
            res += _img(i, j) * k;
            normFactor += k;
        }
    }
    
    return uchar(max(min(res / normFactor, 254.), 0.));
}

    
