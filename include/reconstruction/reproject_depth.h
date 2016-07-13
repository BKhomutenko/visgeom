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
Takes the depthmap from the first image, reconstructs the pointcloud
in the frame of the second image, then reprojects into second image
frame. Reconstructs the cloud at the image points specified by the 
first pointcloud, and sends this back to the first image. This 
pointcloud is transformed into the frame of the first image, and each
point in this cloud is projected onto the line of it's original point.
This new depthmap is the reprojected depthmap.
*/

//TODO: Add image to explain the idea

#pragma once

#include "geometry/geometry.h"
#include "reconstruction/depth_map.h"

class ReprojectedDepth
{
public:
	ReprojectedDepth() {}

	bool reprojectionBetweenFrames(const DepthMap& dMap1, const DepthMap& dMap2,
	        const Transformation<double> T12, DepthMap& output);
private:
	//null
};
