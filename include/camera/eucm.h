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
    bool operator() (const T* params, const T* src, T* dst) const
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
        
        // Check that the point is in the upper hemisphere in case of ellipsoid
        if (alpha > T(0.5))
        {
            const T zn = z / denom; 
            const T C = (alpha - T(1.)) / (alpha + alpha - T(1.));
            if (zn < C) return false;
        }
        
        // Project the point to the mu plane
        const T xn = x / denom;
        const T yn = y / denom;
        // Compute image point
        dst[0] = fu * xn + u0;
        dst[1] = fv * yn + v0;
        return true;  
    } 
      
    static const int INTRINSIC_COUNT = 6;
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
    
    virtual bool projectionJacobian(const Vector3d & src, double * dudx, double * dvdx) const
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

        double rho = sqrt(z * z + beta * (x * x + y * y));
        double gamma = 1. - alpha;
        double eta = alpha * rho + gamma * z;
        
        bool isProjected = true;
        if (eta < 1e-3) isProjected = false;
        else if (alpha > 0.5)
        {
            const double zn = z / eta; 
            const double C = (alpha - 1.) / (alpha + alpha - 1.);
            if (zn < C) isProjected = false;
        }
        
        if (not isProjected)
        {
            dudx[0] = 0;
            dudx[1] = 0;
            dudx[2] = 0;
            dvdx[0] = 0;  
            dvdx[1] = 0;
            dvdx[2] = 0;
            return false;
         
        }
        
        double k = 1. / eta / eta;
        double abrho = alpha * beta / rho;
        double Jxy = k * abrho * x * y;
        double Jz = k * (gamma + alpha * z / rho);
        double Jx = gamma * z + alpha * rho;
        dudx[0] = fu * k * (Jx - abrho * x * x);
        dudx[1] = -fu * Jxy;
        dudx[2] = -fu * x * Jz;
        dvdx[0] = -fv * Jxy;    
        dvdx[1] = fv * k * (Jx - abrho * y * y);
        dvdx[2] = -fv * y * Jz;

        return true;

    }
    
    virtual bool intrinsicJacobian(const Vector3d & src,
                    double * dudalpha, double * dvdalpha) const
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
        
        double x2y2 = x * x + y * y;
        double rho2 = z * z + beta * (x2y2);
        double rho = sqrt(rho2);
        double gamma = 1. - alpha;
        double eta = alpha * rho + gamma * z;
        
        bool isProjected = true;
        if (eta < 1e-3) isProjected = false;
        else if (alpha > 0.5)
        {
            const double zn = z / eta; 
            const double C = (alpha - 1.) / (alpha + alpha - 1.);
            if (zn < C) isProjected = false;
        }
        
        if (not isProjected)
        {
            for (int i = 0; i < 6; i++)
            {
                dudalpha[i] = 0;
                dvdalpha[i] = 0;
            }
            return false;
        }
        
        double eta2 = eta * eta;
        
        dudalpha[0] = -fu * x * (rho - z) / eta2;
        dudalpha[1] = -fu * x * alpha * x2y2 / (2 * eta2 * rho);
        dudalpha[2] = x / eta;
        dudalpha[3] = 0;
        dudalpha[4] = 1;
        dudalpha[5] = 0;
        
        dvdalpha[0] = -fv * y * (rho - z) / eta2;
        dvdalpha[1] = -fv * y * alpha * x2y2 / (2 * eta2 * rho);
        dvdalpha[2] = 0;
        dvdalpha[3] = y / eta;
        dvdalpha[4] = 0;
        dvdalpha[5] = 1;
        
        return true;

    }
    
    virtual double upperBound(int idx) const
    {
        switch (idx)
        {
        case 0:     return 1;       //alpha
        case 1:     return 10;     //beta
        default:    return 1e5;     //the rest
        }
    }
    
    virtual double lowerBound(int idx) const
    {
        switch (idx)
        {
        case 0:     return 0;       //alpha
        case 1:     return 0.1;     //beta
        default:    return 1;     //the rest
        }
    }
    
    virtual EnhancedCamera * clone() const { return new EnhancedCamera(width, height, params.data()); }
    
    virtual ~EnhancedCamera() {}
};

