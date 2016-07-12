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
    double operator() (int x, int y) const
    {
        return (kuu*x + kuv*y + ku)*x 
                + (kvv*y + kv)*y 
                + k1;
    }
    
    double gradx(int x, int y) const
    {
        return 2*kuu*x + kuv*y + ku;
    }
    
    double grady(int x, int y) const
    {
        return kuv*x + 2*kvv*y + + kv;
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
    
    CurveRasterizer(Vector2<T> pt, Vector2<T> epipole, const Surface & surf) :
            x(pt[0]), y(pt[1]), surf(surf)
    {
        fx = surf.gradx(x, y);
        fy = surf.grady(x, y);
        delta = surf(x, y);
        if (fx*(epipole[1] - y) - fy*(epipole[0] - x) > 0) eps = 1;
        else eps = -1;
    }
    
    void moveX(int dx)
    {
        if (dx == 0) return;
        x += dx;
        double fx2 = surf.gradx(x, y);
        delta += 0.5*dx*(fx + fx2);
        fx = fx2;
        fy = surf.grady(x, y); 
    }
    
    void moveY(int dy)
    {
        if (dy == 0) return;
        y += dy;
        double fy2 = surf.grady(x, y);
        delta += 0.5*dy*(fy + fy2);
        fy = fy2;
        fx = surf.gradx(x, y); 
    }
    
    void step()
    {
        if (abs(fx) > abs(fy))  // go along y
        {
            moveY(eps*sign(fx));
            moveX(-round(delta/fx));
        }   
        else  // go along x
        {
            moveX(-eps*sign(fy));
            moveY(-round(delta/fy));
        }
    }
    
    void unstep()
    {
        if (abs(fx) > abs(fy))  // go along y
        {
            moveY(-eps*sign(fx));
            moveX(-round(delta/fx));
        }   
        else  // go along x
        {
            moveX(eps*sign(fy));
            moveY(-round(delta/fy));
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
