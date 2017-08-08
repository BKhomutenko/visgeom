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
Scale parameters. To be used as a parent class wherever scale conversion is needed
NOTE:
(u, v) is an image point 
(x, y) is a depth map point 
*/

#include "reconstruction/scale_parameters.h"


ScaleParameters::ScaleParameters(const ptree & params)
{
    for (auto & item : params)
    {
        const string & pname = item.first;
        if (pname == "scale")       scale = item.second.get_value<int>();
        else if (pname == "u0")     u0 = item.second.get_value<int>();
        else if (pname == "v0")     v0 = item.second.get_value<int>();
        else if (pname == "uMax")   uMax = item.second.get_value<int>();
        else if (pname == "vMax")   vMax = item.second.get_value<int>();
        else if (pname == "xMax")   xMax = item.second.get_value<int>();
        else if (pname == "yMax")   yMax = item.second.get_value<int>();
        else if (pname == "equal_margins" and item.second.get_value<bool>()) setEqualMargin();
    }
//    cout << "u0  " << u0;
//    cout << "v0  " << v0;
//    cout << "xMax  " << xMax;
//    cout << "yMax  " << yMax;
}

void ScaleParameters::setEqualMargin()
{
    xMax = (uMax - 2*u0) / scale + 1;
    yMax = (vMax - 2*v0) / scale + 1;
    if (xMax < 1 or yMax < 1)
    {
        throw std::runtime_error("The scaled image is empty: xMax < 1 or yMax < 1");
    }
}

void ScaleParameters::setXMargin(int val)
{
    xMax = (uMax - u0 - val) / scale + 1;
    if (xMax < 1)
    {
        throw std::runtime_error("The scaled image is empty: xMax < 1");
    }
}

void ScaleParameters::setYMargin(int val)
{
    yMax = (vMax - v0 - val) / scale + 1;
    if (yMax < 1)
    {
        throw std::runtime_error("The scaled image is empty: yMax < 1");
    }
}

// from image to small disparity coordiante transform
int ScaleParameters::xConv(double u) const { return round((u - u0) / scale); }
int ScaleParameters::yConv(double v) const { return round((v - v0) / scale); }

// from small disparity to image coordiante transform    
int ScaleParameters::uConv(int x) const { return x * scale + u0; }
int ScaleParameters::vConv(int y) const { return y * scale + v0; }

bool ScaleParameters::operator == (const ScaleParameters & other) const
{
    return scale == other.scale and u0 == other.u0 and v0 == other.v0 and
            uMax == other.uMax and vMax == other.vMax and 
            xMax == other.xMax and yMax == other.yMax;
}
