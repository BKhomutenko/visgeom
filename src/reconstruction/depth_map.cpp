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
Depth container
*/

#include "reconstruction/depth_map.h"

#include "std.h"
#include "eigen.h"

// nearest neighbor interpolation
double DepthMap::nearest(double u, double v)
{
    return at(x(u), y(v));
}

// nearest neighbor interpolation
double DepthMap::nearest(Vector2d pt)
{
    return at(x(pt[0]), y(pt[1]));
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
double DepthMap::u(int x)
{
    return(x + 0.5)*scale + u0;
}
double DepthMap::v(int y)
{
    return(y + 0.5)*scale + v0;
}

// image coordinates of the block corner
int DepthMap::uc(int x)
{
    return x*scale + u0;
}
int DepthMap::vc(int y)
{
    return y*scale + v0;
}

// depth coordinates of image points
int DepthMap::x(double u)
{
    return floor((u - u0) / scale);
}

int DepthMap::y(double v)
{
    return floor((v - v0) / scale);
}

void DepthMap::reconstruct(Vector3dVec & result)
{
    Vector2dVec pointVec;
    pointVec.reserve(width * height);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            pointVec.emplace_back(x, y);
        }
    }
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < valVec.size(); i++)
    {
        result[i] = result[i].normalized()*valVec[i];
    }
}

void DepthMap::reconstruct(const vector<int> & indexVec, Vector3dVec & result)
{
    Vector2dVec pointVec;
    pointVec.reserve(indexVec.size());
    for (int index : indexVec) // TODO make safe
    {
        int x = index % width;
        int y = index / width;
        pointVec.emplace_back(x, y);
    }
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < indexVec.size(); i++)
    {
        result[i] = result[i].normalized()*valVec[indexVec[i]];
    }
}

void DepthMap::reconstruct(const Vector2dVec & pointVec, Vector3dVec & result)
{
    cameraPtr->reconstructPointCloud(pointVec, result);
    for (int i = 0; i < valVec.size(); i++)
    {
        result[i] = result[i].normalized()*nearest(pointVec[i]);
    }
}

