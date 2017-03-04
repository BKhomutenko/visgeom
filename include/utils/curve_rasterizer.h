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
A template-base coordinate transformation implementation
*/

#pragma once

#include "std.h"
#include "eigen.h"

// Curves in images

// Polynomial with two variables of degree 2
// Axx + Bxy + Cyy + Dx + Ey + F
struct Polynomial2
{
    Polynomial2() {}
//    Polynomial2(const std::array<double, 6> & coeffs) : coeffArr(coeffs) {}
//    std::array<double, 6> coeffArr;
    double kuu, kuv, kvv, ku, kv, k1;
    double operator() (int u, int v) const
    {
        return (kuu*u + kuv*v + ku)*u 
                + (kvv*v + kv)*v 
                + k1;
    }
    
    double gradu(int u, int v) const
    {
        return 2*kuu*u + kuv*v + ku;
    }
    
    double gradv(int u, int v) const
    {
        return kuv*u + 2*kvv*v + + kv;
    }
    
};

inline int sign(double x)
{
    return 2*int(x > 0) - 1;
}

/*
//TODO int or double?
template<typename T, typename Surface>
struct CurveRasterizer
{
    double delta;
    double fx, fy;
    int eps;
    T x, y;
    Surface surf; 
    
    CurveRasterizer(T x, T y, T ex, T ey, const Surface & surf) :
            x(x), y(y), surf(surf)
    {
        fx = surf.gradx(x, y);
        fy = surf.grady(x, y);
        delta = surf(x, y);
        if (fx*(ey - y) - fy*(ex - x) > 0) eps = 1;
        else eps = -1;
    }
    
//    void moveX(int dx)
//    {
//        if (dx == 0) return;
//        x += dx;
//        double fx2 = surf.gradx(x, y);
//        delta += 0.5*dx*(fx + fx2);
//        fx = fx2;
//        fy = surf.grady(x, y); 
//    }
//    
//    void moveY(int dy)
//    {
//        if (dy == 0) return;
//        y += dy;
//        double fy2 = surf.grady(x, y);
//        delta += 0.5*dy*(fy + fy2);
//        fy = fy2;
//        fx = surf.gradx(x, y); 
//    }
    
    void step()
    {
        double fx = surf.gradx(x, y);
        double fy = surf.grady(x, y);
        if (abs(fx) > abs(fy))  // go along y
        {
            y += eps*sign(fx);
            x -= round(surf(x, y)/surf.gradx(x, y));
        }   
        else  // go along x
        {
            x -= eps*sign(fy);
            y -= round(surf(x, y)/surf.grady(x, y));
        }
    }
    
    void unstep()
    {
        double fx = surf.gradx(x, y);
        double fy = surf.grady(x, y);
        if (abs(fx) > abs(fy))  // go along y
        {
            y -= eps*sign(fx);
            x -= round(surf(x, y)/surf.gradx(x, y));
        }   
        else  // go along x
        {
            x += eps*sign(fy);
            y -= round(surf(x, y)/surf.grady(x, y));
        }
    }
    
    void steps(int nsteps)
    {
        if (nsteps > 0)
        {
            for (int i = 0; i < nsteps; i++)
            {
                step();
            }
        }
        else
        {
            for (int i = 0; i > nsteps; i--)
            {
                unstep();
            }
        }
    }
};*/

template<typename T, class Surface>
struct CurveRasterizer
{
    double delta;
    double fu, fv;
    int eps;
    T u, v;
    Surface surf; 
    
    CurveRasterizer(const T u, const T v, const T eu, const T ev, const Surface & surf) :
            u(u), v(v), surf(surf)
    {
        fu = surf.gradu(u, v);
        fv = surf.gradv(u, v);
        delta = surf(u, v);
        if (fu*(ev - v) - fv*(eu - u) > 0) eps = 1;
        else eps = -1;
    }
    
    CurveRasterizer(const Vector2<T> pt, const Vector2<T> epipole, const Surface & surf) :
            u(pt[0]), v(pt[1]), surf(surf)
    {
        fu = surf.gradu(u, v);
        fv = surf.gradv(u, v);
        delta = surf(u, v);
        if (fu*(epipole[1] - v) - fv*(epipole[0] - u) > 0) eps = 1;
        else eps = -1;
    }
    
    void setStep(int step)
    {
        eps *= step;
    }
    
    void moveU(int du)
    {
        if (du == 0) return;
        u += du;
        double fu2 = surf.gradu(u, v);
        delta += 0.5*du*(fu + fu2);
        fu = fu2;
        fv = surf.gradv(u, v); 
    }
    
    void moveV(int dv)
    {
        if (dv == 0) return;
        v += dv;
        double fv2 = surf.gradv(u, v);
        delta += 0.5*dv*(fv + fv2);
        fv = fv2;
        fu = surf.gradu(u, v); 
    }
    
    void step()
    {
        if (abs(fu) > abs(fv))  // go along y
        {
            moveV(eps*sign(fu));
            moveU(-round(delta/fu));
        }   
        else  // go along x
        {
            moveU(-eps*sign(fv));
            moveV(-round(delta/fv));
        }
    }
    
    void unstep()
    {
        if (abs(fu) > abs(fv))  // go along y
        {
            moveV(-eps*sign(fu));
            moveU(-round(delta/fu));
        }   
        else  // go along x
        {
            moveU(eps*sign(fv));
            moveV(-round(delta/fv));
        }
    }
    
    void steps(int nsteps)
    {
        if (nsteps > 0)
        {
            for (int i = 0; i < nsteps; i++)
            {
                step();
            }
        }
        else
        {
            for (int i = 0; i > nsteps; i--)
            {
                unstep();
            }
        }
    }
};
