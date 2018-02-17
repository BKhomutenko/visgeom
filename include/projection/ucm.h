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

#include "projection/generic_camera.h"


//FIXME make a proper unified camera
//so far it is just EUCM with beta = 1

template<typename T> 
struct UnifiedProjector
{
    bool operator() (const T* params, const T* src, T* dst) const
    {
        //TODO check whether projected
        const T & xi = params[0];
        const T & fu = params[1];
        const T & fv = params[2];
        const T & u0 = params[3];
        const T & v0 = params[4];
        
        const T & x = src[0];
        const T & y = src[1];
        const T & z = src[2];
        
        T rho = sqrt(z*z + x*x + y*y);
        T denominv = T(1.) / (z + xi*rho);

        // Project the point to the mu plane
        T xn = x * denominv;
        T yn = y * denominv;
        
        // Compute image point
        dst[0] = fu * xn + u0;
        dst[1] = fv * yn + v0;
        return true;  
    } 
      
    static const int INTRINSIC_COUNT = 5;
};

class UnifiedCamera : public ICamera
{
public:
    using ICamera::params;
    using ICamera::width;
    using ICamera::height;
    UnifiedCamera(int W, int H, const double * const parameters) : ICamera(W, H, 5)
    {  
        ICamera::setParameters(parameters);
    }

    UnifiedCamera(const double * const parameters)  : ICamera(1, 1, 5)
    {  
        ICamera::setParameters(parameters);
    }
    
     /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        //FIXME ???
        const double & xi = params[0];
        const double & fu = params[1];
        const double & fv = params[2];
        const double & u0 = params[3];
        const double & v0 = params[4];
        
        double xn = (src(0) - u0) / fu;
        double yn = (src(1) - v0) / fv;
        
        double u2 = xn * xn + yn * yn;
        
        
        double gamma = sqrt(1. + u2*(1 - xi*xi));
            
        double etanum = -gamma - xi*u2;
        double etadenom = xi*xi*u2 - 1;
        dst << xn, yn, etadenom/(etadenom + xi*etanum);
        
        return true;
    }

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const
    {
        UnifiedProjector<double> projector;
        return projector(params.data(), src.data(), dst.data()); 
    }
    
    virtual bool projectionJacobian(const Vector3d & src, double * dudx, double * dvdx) const
    {
        const double & xi = params[0];
        const double & fu = params[1];
        const double & fv = params[2];
        const double & u0 = params[3];
        const double & v0 = params[4];
        
        const double & x = src[0];
        const double & y = src[1];
        const double & z = src[2];
        
        double xx = x*x;
        double yy = y*y;
        double zz = z*z;
        double rho = sqrt(xx + yy + zz);
        double rhoinv = 1. / rho;
        double deninv = 1./ (xi*rho + z);
        double deninv2 = deninv * deninv;
        
        
        //normalized point
        double xn = x * deninv;
        double yn = y * deninv;
        
        Covector3d dudX; //  normalized point Jacobian dm / dX 
        dudX(0) = (xi*rho + z - xi*xx*rhoinv) * deninv2;
        dudX(1) = -xi * x * y * rhoinv * deninv2;
        dudX(2) = -x * (1 + xi*z*rhoinv) * deninv2;
        
        Covector3d dvdX;
        dvdX(0) = -xi * x * y * rhoinv * deninv2;
        dvdX(1) = (xi*rho + z - xi*yy*rhoinv) * deninv2;
        dvdX(2) = -y * (1 + xi*z*rhoinv) * deninv2;
        
        Map<Covector3d>((double *)dudx) = fu * dudX;
        Map<Covector3d>((double *)dvdx) = fv * dvdX;
        return true;

    }
    
    virtual bool intrinsicJacobian(const Vector3d & src,
                    double * dudalpha, double * dvdalpha) const
    {
        const double & xi = params[0];
        const double & fu = params[1];
        const double & fv = params[2];
        const double & u0 = params[3];
        const double & v0 = params[4];
        
        const double & x = src[0];
        const double & y = src[1];
        const double & z = src[2];
        
        double xx = x*x;
        double yy = y*y;
        double zz = z*z;
        double rho = sqrt(xx + yy + zz);
        double rhoinv = 1. / rho;
        double deninv = 1. / (xi*rho + z);
        double deninv2 = deninv * deninv;
        
        
        //normalized point
        double xn = x * deninv;
        double yn = y * deninv;
        
        //xi
        dudalpha[0] = -fu * xn * deninv * rho;
        //Projection matrix part
        dudalpha[1] = xn;
        dudalpha[2] = 0;
        dudalpha[3] = 1;
        dudalpha[4] = 0;
        
        
        //xi
        dvdalpha[0] = -fv * yn * deninv * rho;
        //Projection matrix part
        dvdalpha[1] = 0;
        dvdalpha[2] = yn;
        dvdalpha[3] = 0;
        dvdalpha[4] = 1;
        return true;

    }
    
    virtual double upperBound(int idx) const
    {
        switch (idx)
        {
        case 0:     return 3;       //xi
        default:    return 1e5;     //the rest
        }
    }
    
    virtual double lowerBound(int idx) const
    {
        switch (idx)
        {
        case 0:     return 0;       //xi
        default:    return 1;     //the rest
        }
    }
    
    virtual double getCenterU() { return params[3]; }
    virtual double getCenterV() { return params[4]; }
    
    virtual UnifiedCamera * clone() const { return new UnifiedCamera(width, height, params.data()); }
    
    virtual ~UnifiedCamera() {}
};

