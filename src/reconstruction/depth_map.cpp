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

//check the limits
bool DepthMap::isValid(int x, int y) const
{
    return (x >= 0 and x < width and y >= 0 and y < height);
}

//FIXME bug, if we write otside the range then read from outside the range, we'll get the written value
// because of foo

// nearest neighbor interpolation
const double & DepthMap::nearest(int u, int v) const
{
    int xd = x(u);
    int yd = y(v);
    if (isValid(xd, yd)) return at(xd, yd);
    else return foo;
}

// nearest neighbor interpolation
const double & DepthMap::nearest(Vector2d pt) const
{
    int xd = x(pt[0]);
    int yd = y(pt[1]);
    if (isValid(xd, yd)) return at(xd, yd);
    else return foo;
}

double & DepthMap::nearest(Vector2d pt)
{
    int xd = x(pt[0]);
    int yd = y(pt[1]);
    if (isValid(xd, yd)) return at(xd, yd);
    else return foo;
}


// nearest neighbor interpolation
const double & DepthMap::nearestSigma(int u, int v) const
{
    int xd = x(u);
    int yd = y(v);
    if (isValid(xd, yd)) return sigma(xd, yd);
    else return foo;
}

// nearest neighbor interpolation
const double & DepthMap::nearestSigma(Vector2d pt) const
{
    int xd = x(pt[0]);
    int yd = y(pt[1]);
    if (isValid(xd, yd)) return sigma(xd, yd);
    else return foo;
}

double & DepthMap::nearestSigma(Vector2d pt)
{
    int xd = x(pt[0]);
    int yd = y(pt[1]);
    if (isValid(xd, yd)) return sigma(xd, yd);
    else return foo;
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

// to access the elements directly
double & DepthMap::at(int idx)
{
    return valVec[idx];
}
const double & DepthMap::at(int idx) const
{
    return valVec[idx];
}

// to access the uncertainty directly
double & DepthMap::sigma(int x, int y)
{
    return sigmaVec[x + y*width];
}
const double & DepthMap::sigma(int x, int y) const
{
    return sigmaVec[x + y*width];
}

// to access the uncertainty directly
double & DepthMap::sigma(int idx)
{
    return sigmaVec[idx];
}
const double & DepthMap::sigma(int idx) const
{
    return sigmaVec[idx];
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

void DepthMap::getPointVec(const std::vector<int> idxVec, Vector2dVec & result) const
{
    result.clear();
    result.reserve(idxVec.size());
    for (auto & idx : idxVec)
    {
        result.emplace_back(u(idx % width), v(idx / width));
    }
}

void DepthMap::reconstructUncertainty(vector<int> & idxVec, 
            Vector3dVec & minDistVec, Vector3dVec & maxDistVec) const
{
    minDistVec.clear();
    maxDistVec.clear();
    idxVec.clear();
    Vector2dVec pointBrutVec;
    vector<double> minVec;
    vector<double> maxVec;
    vector<int> idxBrutVec;
    for (int i = 0; i < valVec.size(); i++)
    {
        double d = valVec[i];
        if (d >= MIN_DEPTH)
        {
            double s = sigmaVec[i];
            // take d +- 2*sigma
            minVec.push_back(max(MIN_DEPTH, d - 2*s));
            maxVec.push_back(d + 2*s);
            idxBrutVec.push_back(i);
        }
    }
    
    getPointVec(idxBrutVec, pointBrutVec);
    
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(pointBrutVec, reconstBrutVec, maskVec);
    
    for (int i = 0; i < reconstBrutVec.size(); i++)
    {
        if (maskVec[i])
        {
            Vector3d X = reconstBrutVec[i].normalized();
            minDistVec.push_back(X*minVec[i]);
            maxDistVec.push_back(X*maxVec[i]);
            idxVec.push_back(idxBrutVec[i]);        
        }
    }
}

void DepthMap::reconstruct(vector<int> & idxVec, Vector3dVec & result) const
{
    result.clear();
    idxVec.clear();
    Vector2dVec pointBrutVec;
    vector<double> depthVec;
    vector<int> idxBrutVec;
    for (int i = 0; i < valVec.size(); i++)
    {
        double d = valVec[i];
        if (d >= MIN_DEPTH)
        {
            depthVec.push_back(d);
            idxBrutVec.push_back(i);
        }
    }
    getPointVec(idxBrutVec, pointBrutVec);
    
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(pointBrutVec, reconstBrutVec, maskVec);
    
    for (int i = 0; i < reconstBrutVec.size(); i++)
    {
        if (maskVec[i])
        {
            Vector3d X = reconstBrutVec[i].normalized();
            result.push_back(X*depthVec[i]);
            idxVec.push_back(idxBrutVec[i]);        
        }
    }
}

void DepthMap::reconstruct(const Vector2dVec & queryPointVec,
        vector<int> & idxVec, Vector3dVec & result) const
{
    result.clear();
    idxVec.clear();
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(queryPointVec, reconstBrutVec, maskVec);
    for (int i = 0; i < queryPointVec.size(); i++)
    {
        if (maskVec[i])
        {
            double d = nearest(queryPointVec[i]);
            if (d < MIN_DEPTH) continue;
            Vector3d X = reconstBrutVec[i].normalized();
            result.push_back(X*d);
            idxVec.push_back(i);        
        }
    }
}

void DepthMap::reconstruct(const Vector2dVec & queryPointVec, Vector3dVec & result) const
{
    vector<int> foo;
    reconstruct(queryPointVec, foo, result);
}

void DepthMap::project(const Vector3dVec & pointVec, Vector2dVec & result) const
{
    cameraPtr->projectPointCloud(pointVec, result);
}

//TODO do not reconstruct all the points but a selected subset
// to avoid reconstruction of points with bad disparity
void DepthReprojector::wrapDepth(const DepthMap& dMap1, const DepthMap& dMap2,
        const Transformation<double> T12, DepthMap& output)
{
	//Step 1 : Get point-cloud of first camera in first frame
	vector<int> idx0Vec;
	Vector3dVec cloud11;
	dMap1.reconstruct(idx0Vec, cloud11);

	//Step 2 : Transform above into second frame
	Vector3dVec cloud12;
	T12.inverseTransform(cloud11, cloud12);

	//Step 3 : Reproject points into second camera
	Vector2dVec point12Vec;
	dMap2.project(cloud12, point12Vec);

	//Step 4 : For reprojected points, reconstruct point-cloud of second camera in second frame
	Vector3dVec cloud22;
	vector<int> idx1Vec;
	dMap2.reconstruct(point12Vec, idx1Vec, cloud22);

	//Step 5 : Transform above into first frame
	Vector3dVec cloud21;
	T12.transform(cloud22, cloud21);

	//Step 6 : Project above points along corresponding depth vectors
    output = dMap1;
    output.setTo(0, 1);
	for(int i = 0; i < idx1Vec.size(); ++i)
	{
	    int idx1 = idx1Vec[i];
	    int idx0 = idx0Vec[idx1];
		const Vector3d & X2 = cloud21[i];
		const Vector3d & X1 = cloud11[idx1];
		//Calculate dot-product to get the distance as the projection along the line
		output.at(idx0) =  X2.dot(X1.normalized());
		output.sigma(idx0) = dMap2.nearestSigma(point12Vec[idx1]);
	}
}

