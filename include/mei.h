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
The Unified Camera Model with distortions
*/

#pragma once

#include <Eigen/Eigen>

#include "camera.h"


template<typename T> 
struct MeiProjector
{
    bool operator () (const T* params, const T* src, T* dst)
    {
        const T & xi = params[0];
        const T & k1 = params[1];
        const T & k2 = params[2];
        const T & k3 = params[3];
        const T & k4 = params[4];
        const T & k5 = params[5];
        const T & fu = params[6];
        const T & fv = params[7];
        const T & u0 = params[8];
        const T & v0 = params[9];
        
        const T & x = src[0];
        const T & y = src[1];
        const T & z = src[2];
        
        T rho = sqrt(z*z + x*x + y*y);
        T denom = z + xi*rho;

        // Project the point to the mu plane
        T xn = x / denom;
        T yn = y / denom;

        // Apply the distortion
        T xx = xn*xn, xy = xn*yn, yy = yn*yn;
        T r2 = xx + yy;
        T D = T(1.) + k1*r2 + k2*r2*r2 + k3*r2*r2*r2;
        
        T deltax = T(2.)*k4*xy + k5*(r2 + T(2.)*xx);
        T deltay = T(2.)*k5*xy + k4*(r2 + T(2.)*yy);
        
        // Compute image point
        dst[0] = fu * (xn*D + deltax) + u0;
        dst[1] = fv * (yn*D + deltay) + v0;
        return true;  
    } 
};

class MeiCamera : public ICamera
{
public:
    using ICamera::params;
    using ICamera::width;
    using ICamera::height;
    MeiCamera(int W, int H, const double * const parameters) : ICamera(W, H, 6)
    {  
        ICamera::setParameters(parameters);
    }

    MeiCamera(const double * const parameters)  : ICamera(1, 1, 6)
    {  
        ICamera::setParameters(parameters);
    }
    
     /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        return false;
    }

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const
    {
        MeiProjector<double> projector;
        return projector(params.data(), src.data(), dst.data());
    }
    
    virtual bool projectionJacobian(const Vector3d & src, Eigen::Matrix<double, 2, 3> & Jac) const
    {
        return false;

    }
    
    
    virtual MeiCamera * clone() const { return new MeiCamera(width, height, params.data()); }
    
    virtual ~MeiCamera() {}
};

