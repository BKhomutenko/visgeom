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




#include "reconstruction/reproject_depth.h"

#include "std.h"

#include "geometry/geometry.h"
#include "reconstruction/depth_map.h"

using std::vector;

bool ReprojectedDepth::reprojectionBetweenFrames(const DepthMap& dMap1, const DepthMap& dMap2,
        const Transformation<double> T12, DepthMap& output)
{
	//Step 1 : Get point-cloud of first camera in first frame
	vector<Vector3d> cloud1;
	dMap1.reconstruct(cloud1);

	//Step 2 : Transform above into second frame
	vector<Vector3d> cloud12;
	T12.transform(cloud1, cloud12);

	//Step 3 : Reproject points into second camera
	vector<Vector2d> reproj;
	dMap2.project(cloud12, reproj);

	//Step 4 : For reprojected points, reconstruct point-cloud of second camera in second frame
	vector<Vector3d> cloud2_filtered;
	dMap2.reconstruct(reproj, cloud2_filtered);

	//Step 5 : Transform above into first frame
	vector<Vector3d> cloud21;
	T12.inverseTransform(cloud2_filtered, cloud21);

	//Step 6 : Project above points along corresponding depth vectors
	vector<double> reconstdist; // Holds the reconstructed distances
	vector<double> reconstsigma; //Holds the reconstructed sigmas
	
	for(int i=0; i<cloud21.size(); ++i)
	{
		const Vector3d point2 = cloud21[i];
		const Vector3d point1 = cloud1[i];
		//Calculate dot-product to get the distance as the projection along the line
		const double lineproj = (point1.normalized()).dot(point2);
		reconstdist.push_back(lineproj);
		const double sigma = dMap2.sigma( dMap2.x(reproj[i][0]), dMap2.y(reproj[i][1]) );
		reconstsigma.push_back( sigma ); // Assume the uncertainty remains same even after transform
	}

	output = dMap1; // Clone the input headers (also clones vectors)
	const int width=dMap1.getWidth(), height=dMap1.getHeight();
	for (int y=0; y<height; ++y)
	{
		for(int x=0; x<width; ++x)
		{
			output.at(x,y) = reconstdist[y*width + x];
			output.sigma(x,y) = reconstsigma[y*width + x];
		}
	}

	return true;
}
