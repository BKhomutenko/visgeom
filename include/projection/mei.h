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

#include "eigen.h"

#include "projection/generic_camera.h"


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
        T denominv = T(1.) / (z + xi*rho);

        // Project the point to the mu plane
        T xn = x * denominv;
        T yn = y * denominv;

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
    
    static const int INTRINSIC_COUNT = 10;
};

class MeiCamera : public ICamera
{
public:
    using ICamera::params;
    using ICamera::width;
    using ICamera::height;
    MeiCamera(int W, int H, const double * const parameters) : ICamera(W, H, 10)
    {  
        ICamera::setParameters(parameters);
    }

    MeiCamera(const double * const parameters)  : ICamera(1, 1, 10)
    {  
        ICamera::setParameters(parameters);
    }
    
     /// takes raw image points and apply undistortion model to them
    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        //FIXME ???
        const double & xi = params[0];
        const double & fu = params[6];
        const double & fv = params[7];
        const double & u0 = params[8];
        const double & v0 = params[9];
        
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
        MeiProjector<double> projector;
        return projector(params.data(), src.data(), dst.data());
    }
    
    virtual bool projectionJacobian(const Vector3d & src, double * dudx, double * dvdx) const
    {
        const double & xi = params[0];
        const double & k1 = params[1];
        const double & k2 = params[2];
        const double & k3 = params[3];
        const double & k4 = params[4];
        const double & k5 = params[5];
        const double & fu = params[6];
        const double & fv = params[7];
        const double & u0 = params[8];
        const double & v0 = params[9];
        
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
        
        Matrix23drm jac_m; //  normalized point Jacobian dm / dX 
        jac_m(0, 0) = (xi*rho + z - xi*xx*rhoinv) * deninv2;
        jac_m(0, 1) = -xi * x * y * rhoinv * deninv2;
        jac_m(0, 2) = -x * (1 + xi*z*rhoinv) * deninv2;
        
        jac_m(1, 0) = -xi * x * y * rhoinv * deninv2;
        jac_m(1, 1) = (xi*rho + z - xi*yy*rhoinv) * deninv2;
        jac_m(1, 2) = -y * (1 + xi*z*rhoinv) * deninv2;
        
        //distorted point
        double xxn = xn * xn;
        double yyn = yn * yn;
        double xyn = xn * yn;
        double r2 = yn * yn + xn * xn;
        double D = 1. + k1*r2 + k2*r2*r2 + k3*r2*r2*r2;
        
//        double deltax = 2.*k4*xyn + k5*(r2 + 2.*xxn);
//        double deltay = 2.*k5*xyn + k4*(r2 + 2.*yyn);
//        
//        double xd = xn * D + deltax;
//        double yd = yn * D + deltay;

        double dDdr2 = k1 + 2*k2*r2 + 3*k3*r2*r2; // dD/d(r^2)
        // d(r^2)/dx = 2x
        // d(r^2)/dy = 2y
        Covector2d dudmn; //  distortion Jacobian dm_d / dm
        dudmn(0) = D + 2*xxn*dDdr2 + 2*k4*yn + 6*k5*xn;
        dudmn(1) = 2*xyn*dDdr2 + 2*k4*xn + 2*k5*yn;
        
        Covector2d dvdmn; //  distortion Jacobian dm_d / dm
        dvdmn(0) = 2*xyn*dDdr2 + 2*k5*yn + 2*k4*xn;
        dvdmn(1) = D + 2*yyn*dDdr2 + 2*k5*xn + 6*k4*yn;
        
        //apply the projection matrix
        dudmn *= fu;
        dvdmn *= fv;
        
        Map<Covector3d>((double *)dudx) = dudmn * jac_m;
        Map<Covector3d>((double *)dvdx) = dvdmn * jac_m;
        return true;
    }
    
    virtual bool intrinsicJacobian(const Vector3d & src,
                    double * dudalpha, double * dvdalpha) const
    {
        const double & xi = params[0];
        const double & k1 = params[1];
        const double & k2 = params[2];
        const double & k3 = params[3];
        const double & k4 = params[4];
        const double & k5 = params[5];
        const double & fu = params[6];
        const double & fv = params[7];
        const double & u0 = params[8];
        const double & v0 = params[9];
        
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
        
        Matrix23drm jac_m; //  normalized point Jacobian dm / dX 
        jac_m(0, 0) = (xi*rho + z - xi*xx*rhoinv) * deninv2;
        jac_m(0, 1) = -xi * x * y * rhoinv * deninv2;
        jac_m(0, 2) = -x * (1 + xi*z*rhoinv) * deninv2;
        
        jac_m(1, 0) = -xi * x * y * rhoinv * deninv2;
        jac_m(1, 1) = (xi*rho + z - xi*yy*rhoinv) * deninv2;
        jac_m(1, 2) = -y * (1 + xi*z*rhoinv) * deninv2;
        
        //distorted point
        double xxn = xn * xn;
        double yyn = yn * yn;
        double xyn = xn * yn;
        double r2 = yn*yn + xn*xn;
        double D = 1. + k1*r2 + k2*r2*r2 + k3*r2*r2*r2;
        double dDdr2 = k1 + 2*k2*r2 + 3*k3*r2*r2; // dD/d(r^2)
        
        double deltax = 2.*k4*xyn + k5*(r2 + 2.*xxn);
        double deltay = 2.*k5*xyn + k4*(r2 + 2.*yyn);
        
        double xd = xn * D + deltax;
        double yd = yn * D + deltay;
        
        //Needed for dp/dxi
        Covector2d dudmn; //  distortion Jacobian dm_d / dm
        dudmn(0) = D + 2*xxn*dDdr2 + 2*k4*yn + 6*k5*xn;
        dudmn(1) = 2*xyn*dDdr2 + 2*k4*xn + 2*k5*yn;
        
        Covector2d dvdmn; //  distortion Jacobian dm_d / dm
        dvdmn(0) = 2*xyn*dDdr2 + 2*k5*yn + 2*k4*xn;
        dvdmn(1) = D + 2*yyn*dDdr2 + 2*k5*xn + 6*k4*yn;
        
        //apply the projection matrix
        dudmn *= fu;
        dvdmn *= fv;
        
        double dxndxi = -xn * deninv * rho;
        double dyndxi = -yn * deninv * rho;
        
        //xi
        dudalpha[0] = dudmn(0) * dxndxi + dudmn(1) * dyndxi;
        //Radial distortion
        dudalpha[1] = fu * xn * r2;
        dudalpha[2] = fu * xn * r2 * r2;
        dudalpha[3] = fu * xn * r2 * r2 * r2;
        //Tangential distortion
        dudalpha[4] = 2.* fu * xyn;
        dudalpha[5] = fu * (r2 + 2.*xxn);
        //Projection matrix part
        dudalpha[6] = xd;
        dudalpha[7] = 0;
        dudalpha[8] = 1;
        dudalpha[9] = 0;
        
        
        //xi
        dvdalpha[0] = dvdmn(0) * dxndxi + dvdmn(1) * dyndxi;
        //Radial distortion
        dvdalpha[1] = fv * yn * r2;
        dvdalpha[2] = fv * yn * r2 * r2;
        dvdalpha[3] = fv * yn * r2 * r2 * r2;
        //Tangential distortion
        dvdalpha[4] = fv * (r2 + 2.*yyn);
        dvdalpha[5] = 2. * fv * xyn;
        //Projection matrix part
        dvdalpha[6] = 0;
        dvdalpha[7] = yd;
        dvdalpha[8] = 0;
        dvdalpha[9] = 1;
        return true;
    }
    
    virtual MeiCamera * clone() const { return new MeiCamera(width, height, params.data()); }
    
    virtual ~MeiCamera() {}
};

