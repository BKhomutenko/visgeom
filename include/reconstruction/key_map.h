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
Structure for maintaining the map as a multi-keyframe approach
*/

#pragma once

#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "geometry/geometry.h"

#include "reconstruction/depth_map.h"

class KeyMap
{
public:
	KeyMap();
	~KeyMap();

	//Functions required:
	// Add depthmap to map
	// Compare depthmaps
	// Find closest depthmap
	// Transform depthmaps between poses
	// etc.

	// Add a depthmap into the keyframe map structure at specified pose
	void addDepthMap(const DepthMap & dmap, const Transformation<double> & pose);

	// Compare a depthmap with depthmap at specified index in the map structure
	double compare(const DepthMap & dmap, const Transformation<double> & pose, const int index) const;

	// Find the depthmap (by index) that's closest in description to the given depthmap
	int findClosestDepthmap(const DepthMap & dmap, const Transformation<double> & pose) const;

	// Data access functions
	const DepthMap & getDepthMap(const int index) const { return depthmapVec[index]; }
	DepthMap & getDepthMap(const int index) { return depthmapVec[index]; }
	const Transformation<double> & getPose(const int index) const { return basePoseVec[index]; }

private:
	std::vector< DepthMap > depthmapVec;
	std::vector< Transformation<double> > basePoseVec;
}