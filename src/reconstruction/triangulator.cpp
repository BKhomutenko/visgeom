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

#include "reconstruction/triangulator.h"

#include "eigen.h"
#include "geometry/geometry.h"
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
void Triangulator::compute(const Vector3d p, Vector3d q,
        double * res1, double * res2,
        double * jac1, double * jac2) const
{
    assert(res1 != NULL or res2 != NULL);
    if (jac1 != NULL)
    {
        assert(res1 != NULL);
    }
    if (jac2 != NULL)
    {
        assert(res2 != NULL);
    }
    q = R * q;
    double pq = p.dot(q);
    double pp = p.dot(p);
    double qq = q.dot(q);
    double tp = t.dot(p);
    double tq = t.dot(q);
    double ppqq = pp * qq;
    double delta = -pp * qq + pq * pq;
    if (abs(delta) < ppqq * 1e-6) // max 100m for 10cm base
    {
        if (res1 != NULL) *res1 = 0.;
        if (res2 != NULL) *res2 = 0.;
        if (jac1 != NULL)
        {
            fill(jac1, jac1 + 6, 0.); 
        }
        if (jac2 != NULL)
        {
            fill(jac2, jac2 + 6, 0.); 
        }
    }
    else
    {
        double deltaInv = 1. / delta;
        double delta1, delta2;
        if (res1 != NULL) 
        {
            delta1 = (-tp * qq + tq * pq);
            *res1 = delta1 * deltaInv;
            
        }
        if (res2 != NULL)
        {
            delta2 = (pp * tq - pq * tp);
            *res2 = delta2 * deltaInv;
        }
        
        if (jac1 != NULL or jac2 != NULL)
        {
            Matrix3d qSkew = hat(q);
            Vector3d jacDeltaOmega = 2 * pq * qSkew * p;
            if (jac1 != NULL)
            {
                Vector3d jacDelta1V = q*pq - p*qq;
                Vector3d jacDelta1Omega = qSkew * (pq * t + tq * p);
                Map<Vector3d> jacLambda1V(jac1);
                Map<Vector3d> jacLambda1Omega(jac1 + 3);
                jacLambda1V = deltaInv * jacDelta1V;
                jacLambda1Omega = deltaInv * (jacDelta1Omega - *res1 * jacDeltaOmega);
            }
            if (jac2 != NULL)
            {
                Vector3d jacDelta2V = q*pp - p*pq;
                Vector3d jacDelta2Omega = qSkew * (pp * t - tp * p);
                Map<Vector3d> jacLambda2V(jac2);
                Map<Vector3d> jacLambda2Omega(jac2 + 3);
                jacLambda2V = deltaInv * jacDelta2V;
                jacLambda2Omega = deltaInv * (jacDelta2Omega - *res2 * jacDeltaOmega);
            }
        }
    }
}

// num / denom if denom is big enough
// linear extrapolation otherwise
double Triangulator::regDiv(double num, double denom) const
{
    if (denom > eps * num)
    {
        return num / denom;
    }
    else if (num == 0) // singularity, very unlikely
    {
        return 2. / eps;
    }
    else
    {
        return 2. / eps - denom / (num * eps * eps);
    }
}

/*
Triangulation is done by solving:

{ t.(l1 p - l2 q - t) = 0
{ (p + q).(l1 p - l2 q - t) = 0

p, q are direction vectors;
t is the translation taken from the transformation

in case if p and q are close to parallel, the reconstruction is regularized

jac is a kinematic jacobian  l_dot = jac  *   (  v  )
                                              (omega)
                                              
Memory must be allocated
*/
void Triangulator::computeRegular(const Vector3d p, Vector3d q,
        double * res1, double * res2,
        double * jac1, double * jac2) const
{
    assert(res1 != NULL or res2 != NULL);
    if (jac1 != NULL)
    {
        assert(res1 != NULL);
    }
    if (jac2 != NULL)
    {
        assert(res2 != NULL);
    }
    q = R * q;
    const Vector3d r = p + q;
    double pq = p.dot(q);
    double pp = p.dot(p);
    double qq = q.dot(q);
    double tp = t.dot(p);
    double tq = t.dot(q);
    double tr = t.dot(r);
    double tt = t.dot(t);
    double rp = r.dot(p);
    double rq = r.dot(q);
    double delta = tp*rq - tq*rp;
    double delta1, delta2;
    if (res1 != NULL) 
    {
        delta1 = tt * rq - tr * tq;
        *res1 = regDiv(delta1, delta); // max 33m

    }
    if (res2 != NULL)
    {
        delta2 = tt * rp - tr * tp;
        *res2 = regDiv(delta2, delta); // max 33m
    }
    
    if (jac1 != NULL or jac2 != NULL)
    {
        double deltaInv = 1. / delta;
        Matrix3d qSkew = hat(q);
        Vector3d jacDeltaV = rq * p - rp * q;
        Vector3d jacDeltaOmega = qSkew * (tp * p - tq * p - rp * t);
        if (jac1 != NULL)
        {
            Vector3d jacDelta1V = 2*rq * t - tq * r - tr * q;
//            Vector3d jacDelta1Omega = qSkew * (tt * p - 2 * tq * t);
            Vector3d jacDelta1Omega = qSkew * (tt * p - tq * t - tr * t);
            Map<Vector3d> jacLambda1V(jac1);
            Map<Vector3d> jacLambda1Omega(jac1 + 3);
            if (delta > eps * delta1)
            {
                jacLambda1V = deltaInv * (jacDelta1V - *res1 * jacDeltaV);
                jacLambda1Omega = deltaInv * (jacDelta1Omega - *res1 * jacDeltaOmega);
            }
            else if (delta1 == 0) //very unlikely
            {
                jacLambda1V = Vector3d(0, 0, 0);
                jacLambda1Omega = Vector3d(0, 0, 0);
            }
            else
            {
                const double coef = 1. / (eps * eps);
                double deltaInv1 = 1. / delta1;
                double k = -coef * deltaInv1;
                double k1 = coef * delta * deltaInv1 * deltaInv1;
                jacLambda1V = k * jacDeltaV + k1 * jacDelta1V;
                jacLambda1Omega = k * jacDeltaOmega + k1 * jacDelta1Omega;
            }
        }
        if (jac2 != NULL)
        {
            Vector3d jacDelta2V = 2*rp * t - tp * r - tr * p;
            Vector3d jacDelta2Omega = qSkew * (tt * p - tp * t);
            Map<Vector3d> jacLambda2V(jac2);
            Map<Vector3d> jacLambda2Omega(jac2 + 3);
            if (delta > eps * delta2)
            {
                jacLambda2V = deltaInv * (jacDelta2V - *res2 * jacDeltaV);
                jacLambda2Omega = deltaInv * (jacDelta2Omega - *res2 * jacDeltaOmega);
            }
            else if (delta2 == 0) //very unlikely
            {
                jacLambda2V = Vector3d(0, 0, 0);
                jacLambda2Omega = Vector3d(0, 0, 0);
            }
            else
            {
                const double coef = -1. / (eps * eps);
                double deltaInv2 = 1. / delta2;
                double k = coef * deltaInv2;
                double k2 = coef * delta * deltaInv2 * deltaInv2;
                jacLambda2V = k * jacDeltaV - k2 * jacDelta2V;
                jacLambda2Omega = k * jacDeltaOmega - k2 * jacDelta2Omega;
            }
        }
    }
}

