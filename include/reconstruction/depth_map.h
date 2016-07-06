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

#include "std.h"
#include "eigen.h"

#include "camera/eucm.h"

class DepthMap
{
public:
    DepthMap(const EnhancedCamera & camera, int w, int h, double u0, double v0, double scale) :
            camera(camera), width(w), height(h), u0(u0), v0(v0), scale(scale),
            valVec(w*h, 1),  sigmaVec(w*h, 100)  {}
    virtual ~DepthMap() {}
    
    // nearest neighbor interpolation
    double nearest (double u, double v);
    double nearest (Vector2d pt);
    
    // to access the elements directly
    double & at (int x, int y);
    const double & at (int x, int y) const;
    
    // image coordinates of depth points
    double u (int x);
    double v (int y);
    
    // depth coordinates of image points
    double x (int u);
    double y (int v);
    
    void reconstruct(Vector3dVec & result);
    void reconstruct(const Vector2dVec & pointVec, Vector3dVec & result);
    
private:
    EnhancedCamera camera;
    int width;
    int height;
    std::vector<double> valVec;
    std::vector<double> sigmaVec; // uncertainty
    
    double u0, v0; // image coordinates of the [0, 0] point
    double scale; // normally > 1, x = (u - u0) / ration 
};
