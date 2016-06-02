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
The Enhanced Unified Camera Model
*/

#pragma once

#include <Eigen/Eigen>

#include "camera/generic_camera.h"


template<typename T> 
struct EnhancedProjector
{
    bool operator() (const T* params, const T* src, T* dst)
    {
        const T & alpha = params[0];
        const T & beta = params[1];
        const T & fu = params[2];
        const T & fv = params[3];
        const T & u0 = params[4];
        const T & v0 = params[5];
        
        const T & x = src[0];
        const T & y = src[1];
        const T & z = src[2];
        
        T denom = alpha * sqrt(z*z + beta*(x*x + y*y)) + (T(1.) - alpha) * z;

        if (denom < 1e-3) return false;
        // Project the point to the mu plane
        T xn = x / denom;
        T yn = y / denom;

        // Compute image point
        dst[0] = fu * xn + u0;
        dst[1] = fv * yn + v0;
        return true;  
    } 
};

class EnhancedCamera : public ICamera
{
public:
    using ICamera::params;
    using ICamera::width;
    using ICamera::height;
    EnhancedCamera(int W, int H, const double * const parameters) : ICamera(W, H, 6)
    {  
        ICamera::setParameters(parameters);
    }

    EnhancedCamera(const double * const parameters)  : ICamera(1, 1, 6)
    {  
        ICamera::setParameters(parameters);
    }
    
     /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        const double & alpha = params[0];
        const double & beta = params[1];
        const double & fu = params[2];
        const double & fv = params[3];
        const double & u0 = params[4];
        const double & v0 = params[5];
        
        double xn = (src(0) - u0) / fu;
        double yn = (src(1) - v0) / fv;
        
        double u2 = xn * xn + yn * yn;
        double gamma = 1. - alpha;    
        double num = 1. - u2 * alpha * alpha * beta;
        double det = 1 - (alpha - gamma)*beta*u2;
        if (det < 0) return false;
        double denom = gamma + alpha*sqrt(det);
        dst << xn, yn, num/denom;

        return true;
    }

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const
    {
        EnhancedProjector<double> projector;
        return projector(params.data(), src.data(), dst.data()); 
    }
    
    virtual bool projectionJacobian(const Vector3d & src, Eigen::Matrix<double, 2, 3> & Jac) const
    {
        const double & alpha = params[0];
        const double & beta = params[1];
        const double & fu = params[2];
        const double & fv = params[3];
        const double & u0 = params[4];
        const double & v0 = params[5];
        
        const double & x = src(0);
        const double & y = src(1);
        const double & z = src(2);

        double rho = sqrt(z*z + beta*(x*x + y*y));
        double gamma = 1. - alpha;
        double d = alpha * rho + gamma * z;
        double k = 1. / d / d;
        Jac(0,0) = fu * k * (gamma * z + alpha * rho - alpha * beta * x * x / rho);
        Jac(0,1) = -fu * k * alpha * beta * x * y / rho;
        Jac(0,2) = -fu * k * x * (gamma + alpha * z / rho);
        Jac(1,0) = -fv * k * alpha * beta * x * y / rho;    
        Jac(1,1) = fv * k * (gamma * z + alpha * rho - alpha * beta * y * y / rho);
        Jac(1,2) = -fv * k * y * (gamma + alpha * z / rho);

        return true;

    }
    
    
    virtual EnhancedCamera * clone() const { return new EnhancedCamera(width, height, params.data()); }
    
    virtual ~EnhancedCamera() {}
};

