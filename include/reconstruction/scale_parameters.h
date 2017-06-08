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

#include "json.h"


struct ScaleParameters
{
    ScaleParameters() {} //TODO remove, for backward compatibility purposes
    ScaleParameters(const ptree & params);
    // basic parameters
    int scale = 1;
    int u0 = 0, v0 = 0;
    int uMax = 1, vMax = 1; //image size
    int xMax = 1, yMax = 1; //scaled size
    
    void setEqualMargin();
    
    void setXMargin(int val);
    
    void setYMargin(int val);
    
    // from image to small disparity coordiante transform
    int xConv(double u) const;
    int yConv(double v) const;
    
    // from small disparity to image coordiante transform    
    int uConv(int x) const;
    int vConv(int y) const;
    
    bool operator == (const ScaleParameters & other) const;
};

