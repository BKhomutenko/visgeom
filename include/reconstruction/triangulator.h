/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUdouble ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 


#pragma once

#include "eigen.h"
#include "std.h"
#include "geometry/geometry.h"


// num / denom if denom is big enough
// linear extrapolation otherwise
double regDiv(double num, double denom, double EPS);


class Triangulator
{
public:
    Triangulator(double epsilon = 3e-2) : eps(epsilon) {}
    Triangulator(const Transf & transf, double epsilon = 1e-3) :
            R(transf.rotMat()), t(transf.trans()), eps(epsilon) {}
    
    virtual ~Triangulator() {}
    
    void setTransformation(const Transf & transf) 
    {
        R = transf.rotMat();
        t = transf.trans();
    }
    
    void setEpsilon(double newEpsilon) 
    {
        eps = newEpsilon;
    }
    
    double regDiv(double num, double denom) const;
       
    /*

    triangulation is done by solving

    { p (l1 p - l2 q - t) = 0
    { q (l1 p - l2 q - t) = 0

    p, q are direction vectors;
    t is a translation taken from the transformation

    jac is a kinematic jacobian  l_dot = jac  *   (  v  )
                                                  (omega)
    Memory must be allocated
    */
    void compute(const Vector3d p, Vector3d q,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL) const;
    
    /*
    Triangulation is done by solving:

    { t.(l1 p - l2 q - t) = 0
    { (p + q).(l1 p - l2 q - t) = 0

    p, q are direction vectors;
    t is the translation taken from the transformation

    in case if p and q are close to parallel, the reconstruction is regularized
    */     
    void computeRegular(const Vector3d p, Vector3d q,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL) const;
    
    void compute(const Vector3dVec & pVec, const Vector3dVec & qVec,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL) const
    {
        assert(pVec.size() == qVec.size());
        assert(res1 != NULL or res2 != NULL);
        if (jac1 != NULL)
        {
            assert(res1 != NULL);
        }
        if (jac2 != NULL)
        {
            assert(res2 != NULL);
        }
        for (int i = 0; i < pVec.size(); i++)
        {
            compute(pVec[i], qVec[i],
                    res1 != NULL ? res1 + i : NULL,
                    res2 != NULL ? res2 + i : NULL,
                    jac1 != NULL ? jac1 + 6 * i : NULL,
                    jac2 != NULL ? jac2 + 6 * i : NULL);
        }
    }
    
    void computeRegular(const Vector3dVec & pVec, const Vector3dVec & qVec,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL) const
    {
        assert(pVec.size() == qVec.size());
        assert(res1 != NULL or res2 != NULL);
        if (jac1 != NULL)
        {
            assert(res1 != NULL);
        }
        if (jac2 != NULL)
        {
            assert(res2 != NULL);
        }
        for (int i = 0; i < pVec.size(); i++)
        {
            computeRegular(pVec[i], qVec[i],
                    res1 != NULL ? res1 + i : NULL,
                    res2 != NULL ? res2 + i : NULL,
                    jac1 != NULL ? jac1 + 6 * i : NULL,
                    jac2 != NULL ? jac2 + 6 * i : NULL);
        }
    }
    
private:
    Matrix3d R;
    Vector3d t;
    double eps;
};

        
