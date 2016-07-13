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
Depth container. 
NOTE:
(u, v) is an image point 
(x, y) is a depth map point
*/

#include "reconstruction/depth_map.h"

#include "io.h"
#include "std.h"
#include "eigen.h"

// nearest neighbor interpolation
double DepthMap::nearest(int u, int v) const
{
    int xd = x(u);
    int yd = y(v);
    if (isValid(xd, yd)) return at(xd, yd);
    else return 0;
}

bool DepthMap::isValid(int x, int y) const
{
    return (x >= 0 and x < width and y >= 0 and y < height);
}

// nearest neighbor interpolation
double DepthMap::nearest(Vector2d pt) const
{
    int xd = x(pt[0]);
    int yd = y(pt[1]);
    if (isValid(xd, yd)) return at(xd, yd);
    else return 0;
}

// to access the elements directly
double & DepthMap::at(int x, int y)
{
    return valVec[x + y*width];
}
const double & DepthMap::at(int x, int y) const
{
    return valVec[x + y*width];
}

// to access the uncertainty directly
double & DepthMap::sigma(int x, int y)
{
    return valVec[x + y*width];
}
const double & DepthMap::sigma(int x, int y) const
{
    return valVec[x + y*width];
}

// image coordinates of depth points
int DepthMap::u(int x) const
{
    return x * scale + u0;
}
int DepthMap::v(int y) const
{
    return y * scale + v0;
}

// depth coordinates of image points
int DepthMap::x(int u) const
{
    return round(double(u - u0) / scale);
}

int DepthMap::y(int v) const
{
    return round(double(v - v0) / scale);
}

void DepthMap::reconstruct(Vector3dVec & result) const
{
    Vector2dVec pointVec;
    pointVec.reserve(width * height);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            pointVec.emplace_back(u(x), v(y));
        }
    }
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < valVec.size(); i++)
    {
        result[i] = result[i].normalized()*valVec[i];
    }
}

void DepthMap::reconstruct(const vector<int> & indexVec, Vector3dVec & result) const
{
    Vector2dVec pointVec;
    pointVec.reserve(indexVec.size());
    for (int index : indexVec) // TODO make safe
    {
        int x = index % width;
        int y = index / width;
        pointVec.emplace_back(u(x), v(y));
    }
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < indexVec.size(); i++)
    {
        result[i] = result[i].normalized()*valVec[indexVec[i]];
    }
}

void DepthMap::reconstruct(const Vector2dVec & pointVec, Vector3dVec & result) const
{
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < pointVec.size(); i++)
    {
        result[i] = result[i].normalized() * nearest(pointVec[i]);
    }
}

void DepthMap::project(const Vector3dVec & pointVec, Vector2dVec & result) const
{
    cameraPtr->projectPointCloud(pointVec, result);
}

