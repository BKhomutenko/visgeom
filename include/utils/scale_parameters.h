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

#pragma once

struct ScaleParameters
{
    // basic parameters
    int scale = 1;
    int u0 = 0, v0 = 0;
    int uMax, vMax;
    int xMax, yMax;
    
    void setEqualMargin()
    {
        xMax = (uMax - 2*u0) / scale + 1;
        yMax = (vMax - 2*v0) / scale + 1;
        if (xMax < 1 or yMax < 1)
        {
            throw std::runtime_error("The scaled image is empty: xMax < 1 or yMax < 1");
        }
    }
    
    // from image to small disparity coordiante transform
    int x(double u) const { return round((u - u0) / scale); }
    int y(double v) const { return round((v - v0) / scale); }
    
    // from small disparity to image coordiante transform    
    int u(int x) const { return x * scale + u0; }
    int v(int y) const { return y * scale + v0; }
};

