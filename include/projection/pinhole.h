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
The Pinhole Camera Model
*/

#pragma once

#include <Eigen/Eigen>

#include "projection/generic_camera.h"

class Pinhole : public ICamera
{
public:

    Pinhole(double u0, double v0, double f)
    : ICamera(2*u0, 2*v0, 3) 
    {
        params[0] = u0;
        params[1] = v0;
        params[2] = f;
    }
    virtual ~Pinhole() {}

    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & u = src(0);
        const double & v = src(1);
        dst << (u - u0)/f, (v - v0)/f, 1;
        return true;
    }

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & x = src(0);
        const double & y = src(1);
        const double & z = src(2);
        if (z < 1e-2)
        {
            dst << -1, -1;
            return false;
        }
        dst << x * f / z + u0, y * f / z + v0;
        return true;
    }

    //TODO implement the projection and distortion Jacobian
    virtual bool projectionJacobian(const Vector3d & src, Matrix23d & Jac) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & x = src(0);
        const double & y = src(1);
        const double & z = src(2);
        double zz = z * z;
        Jac(0, 0) = f/z;
        Jac(0, 1) = 0;
        Jac(0, 2) = -x * f/ zz;
        Jac(1, 0)= 0;
        Jac(1, 1) = f/z;
        Jac(1, 2) = -y * f/ zz;
    }
    
    virtual Pinhole * clone() const
    {
        return new Pinhole(params[0], params[1], params[2]);
    }
};

