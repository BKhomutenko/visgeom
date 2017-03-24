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
reading out .json files
*/

#pragma once

#include "std.h"
#include "except.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "geometry/geometry.h"

using boost::property_tree::ptree;
using boost::property_tree::read_json;

template<typename T>
Transformation<T> transformFromData(const vector<T> & valVec)
{
    if (valVec.size() == 3) //x, y, theta
    {
        return Transformation<T>(valVec[0], valVec[1], 0, 0, 0, valVec[2]);
    }
    else  if (valVec.size() == 6)   //full 6 dof
    {
        return Transformation<T>(valVec.data());
    }
    
    else  if (valVec.size() == 12)    //homogeneous transformation
    {
        Matrix3<T> R;
        R << valVec[0], valVec[1], valVec[2],
             valVec[4], valVec[5], valVec[6],
             valVec[8], valVec[9], valVec[10];
        
        Vector3<T> t(valVec[3], valVec[7], valVec[11]);
        return Transformation<T>(t, R);
    }
    else
    {
        throw runtime_error("invalid trasformation format. must be 3, 6, or 12 values; " 
            + to_string(valVec.size()) + " are given."); 
    }
}

inline vector<double> readVector(const ptree & node)
{
    vector<double> valVec;
    for (auto & x : node)
    {
        valVec.push_back(x.second.get_value<double>());
    }
    return valVec;
}

inline Transformation<double> readTransform(const ptree & node)
{
    return transformFromData(readVector(node));
}

