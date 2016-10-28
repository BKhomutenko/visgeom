
#include "reconstruction/key_map.h"

#include "io.h"

KeyMap::KeyMap()
{
	// Do nothing
}


KeyMap::~KeyMap()
{
	// Do nothing
}


void KeyMap::addDepthMap(const DepthMap & dmap, const Transformation<double> & pose)
{
	depthmapVec.push_back(dmap);
	basePoseVec.push_back(pose);
}


double KeyMap::compare(const DepthMap & dmap. const Transformation<double> & pose, const int index) const
{
	//TODO complete this function
}


int findClosestDepthmap(const DepthMap & dmap, const Transformation<double> & pose) const
{
	int minindex = -1;
	double minvalue = 1e7;

	for (int i = 0; i < depthmapVec.size(); ++i)
	{
		const double value = compare(dmap, pose, i);
		if (value < minvalue)
		{
			minvalue = value;
			minindex = i
		}
	}

	if (minindex == -1)
	{
		std::cerr << "ERROR : Bad data structure in KeyMap" << std::endl;
		std::exit(1);
	}

	return minindex;
}
