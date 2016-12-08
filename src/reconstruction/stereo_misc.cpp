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

#include "reconstruction/stereo_misc.h"

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
void triangulate(const Transf xi, 
        const Vector3dVec & xVec1, const Vector3dVec & xVec2,
        double * res1, double * res2, double * jac1, double * jac2)
{
    assert(xVec1.size() == xVec2.size());
    assert(res1 != NULL or res2 != NULL);
    if (jac1 != NULL)
    {
        assert(res1 != NULL);
    }
    if (jac2 != NULL)
    {
        assert(res2 != NULL);
    }
    Matrix3d R = xi.rotMat();
    for (int i = 0; i < xVec1.size(); i++)
    {
        const Vector3d & p = xVec1[i];
        const Vector3d q = R * xVec2[i];
        const Vector3d t = xi.trans();
        double pq = p.dot(q);
        double pp = p.dot(p);
        double qq = q.dot(q);
        double tp = t.dot(p);
        double tq = t.dot(q);
        double ppqq = pp * qq;
        double delta = -pp * qq + pq * pq;
        if (abs(delta) < ppqq * 1e-6) // max 100m for 10cm base
        {
            if (res1 != NULL) res1[i] = 0.;
            if (res2 != NULL) res2[i] = 0.;
            if (jac1 != NULL)
            {
                fill(jac1 + i * 6, jac1 + (i + 1) * 6, 0.); 
            }
            if (jac2 != NULL)
            {
                fill(jac2 + i * 6, jac2 + (i + 1) * 6, 0.); 
            }
        }
        else
        {
            double deltaInv = 1. / delta;
            double delta1, delta2;
            if (res1 != NULL) 
            {
                delta1 = (-tp * qq + tq * pq);
                res1[i] = delta1 * deltaInv;
                
            }
            if (res2 != NULL)
            {
                delta2 = (pp * tq - pq * tp);
                res2[i] = delta2 * deltaInv;
            }
            
            if (jac1 != NULL or jac2 != NULL)
            {
                Matrix3d qSkew = hat(q);
                Vector3d jacDeltaOmega = 2 * pq * qSkew * p;
                if (jac1 != NULL)
                {
                    Vector3d jacDelta1V = q*pq - p*qq;
                    Vector3d jacDelta1Omega = qSkew * (pq * t + tq * p);
                    Map<Vector3d> jacLambda1V(jac1 + i * 6);
                    Map<Vector3d> jacLambda1Omega(jac1 + i * 6 + 3);
                    jacLambda1V = deltaInv * jacDelta1V;
                    jacLambda1Omega = deltaInv * (jacDelta1Omega - res1[i] * jacDeltaOmega);
                }
                if (jac2 != NULL)
                {
                    Vector3d jacDelta2V = q*pp - p*pq;
                    Vector3d jacDelta2Omega = qSkew * (pp * t - tp * p);
                    Map<Vector3d> jacLambda2V(jac2 + i * 6);
                    Map<Vector3d> jacLambda2Omega(jac2 + i * 6 + 3);
                    jacLambda2V = deltaInv * jacDelta2V;
                    jacLambda2Omega = deltaInv * (jacDelta2Omega - res2[i] * jacDeltaOmega);
                }
            }
        }
    }
}

// num / denom if denom is big enough
// linear extrapolation otherwise
double regDiv(double num, double denom, double EPS)
{
    if (denom > EPS * num)
    {
        return num / denom;
    }
    else if (num == 0) // singularity, very unlikely
    {
        return 2. / EPS;
    }
    else
    {
        return 2. / EPS - denom / (num * EPS * EPS);
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
void triangulateRegular(const Transf xi, 
        const Vector3dVec & xVec1, const Vector3dVec & xVec2,
        double * res1, double * res2, double * jac1, double * jac2)
{
    const double EPS(1. / 33.); // max 33m
    assert(xVec1.size() == xVec2.size());
    assert(res1 != NULL or res2 != NULL);
    if (jac1 != NULL)
    {
        assert(res1 != NULL);
    }
    if (jac2 != NULL)
    {
        assert(res2 != NULL);
    }
    Matrix3d R = xi.rotMat();
    for (int i = 0; i < xVec1.size(); i++)
    {
        const Vector3d & p = xVec1[i];
        const Vector3d q = R * xVec2[i];
        const Vector3d r = p + q;
        const Vector3d t = xi.trans();
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
            res1[i] = regDiv(delta1, delta, EPS); // max 33m

        }
        if (res2 != NULL)
        {
            delta2 = tt * rp - tr * tp;
            res2[i] = regDiv(delta2, delta, EPS); // max 33m
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
                Vector3d jacDelta1Omega = qSkew * (tt * p - 2 * tq * t);
                Map<Vector3d> jacLambda1V(jac1 + i * 6);
                Map<Vector3d> jacLambda1Omega(jac1 + i * 6 + 3);
                if (delta > EPS * delta1)
                {
                    jacLambda1V = deltaInv * (jacDelta1V - res1[i] * jacDeltaV);
                    jacLambda1Omega = deltaInv * (jacDelta1Omega - res1[i] * jacDeltaOmega);
                }
                else if (delta1 == 0) //very unlikely
                {
                    jacLambda1V = Vector3d(0, 0, 0);
                    jacLambda1Omega = Vector3d(0, 0, 0);
                }
                else
                {
                    const double coef = -1. / (EPS * EPS);
                    double deltaInv1 = 1. / delta1;
                    double k = coef * deltaInv1;
                    double k1 = coef * delta * deltaInv1 * deltaInv1;
                    jacLambda1V = k * jacDeltaV - k1 * jacDelta1V;
                    jacLambda1Omega = k * jacDeltaOmega - k1 * jacDelta1Omega;
                }
            }
            if (jac2 != NULL)
            {
                Vector3d jacDelta2V = 2*rp * t - tp * r - tr * p;
                Vector3d jacDelta2Omega = qSkew * (tt * p - 2 * tp * t);
                Map<Vector3d> jacLambda2V(jac2 + i * 6);
                Map<Vector3d> jacLambda2Omega(jac2 + i * 6 + 3);
                if (delta > EPS * delta2)
                {
                    jacLambda2V = deltaInv * (jacDelta2V - res2[i] * jacDeltaV);
                    jacLambda2Omega = deltaInv * (jacDelta2Omega - res2[i] * jacDeltaOmega);
                }
                else if (delta2 == 0) //very unlikely
                {
                    jacLambda2V = Vector3d(0, 0, 0);
                    jacLambda2Omega = Vector3d(0, 0, 0);
                }
                else
                {
                    const double coef = -1. / (EPS * EPS);
                    double deltaInv2 = 1. / delta2;
                    double k = coef * deltaInv2;
                    double k2 = coef * delta * deltaInv2 * deltaInv2;
                    jacLambda2V = k * jacDeltaV - k2 * jacDelta2V;
                    jacLambda2Omega = k * jacDeltaOmega - k2 * jacDelta2Omega;
                }
            }
        }
    }
}


int computeError(int v, int thMin, int thMax)
{
    return max(0, max(thMin - v, v - thMax));
}

vector<int> compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec, int flawCost)
{
    vector<int> thMinVec(desc.size()), thMaxVec(desc.size());
    
    for (int i = 1; i < desc.size() - 1; i++)
    {
        const int & d = desc[i];
        int d1 = (desc[i] + desc[i - 1]) / 2;
        int d2 = (desc[i] + desc[i + 1]) / 2;
        thMinVec[i] = min(d, min(d1, d2));
        thMaxVec[i] = max(d, max(d1, d2));
    }
    
    if (desc[0] > desc[1])
    {
        thMinVec[0] = (desc[0] + desc[1]) / 2;
        thMaxVec[0] = desc[0]; 
    }
    else
    {
        thMaxVec[0] = (desc[0] + desc[1]) / 2;
        thMinVec[0] = desc[0];
    }
    
    const int & d = desc[desc.size() - 2];
    if (desc.back() > d)
    {
        thMinVec.back() = (desc.back() + d) / 2;
        thMaxVec.back() = desc.back(); 
    }
    else
    {
        thMaxVec.back() = (desc.back() + d) / 2;
        thMinVec.back() = desc.back();
    }
    
    
    const int HALF_LENGTH = desc.size() / 2;
    vector<int> rowA(sampleVec.size()), rowB(sampleVec.size());
    
    //match the first half
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(int(sampleVec[i]) - int(desc[0]));
        rowA[i] = computeError(sampleVec[i], thMinVec[0], thMaxVec[0]);
    }
    for (int i = 1; i <= HALF_LENGTH; i++)
    {
//        rowB[0] = rowA[0] + flawCost + abs(int(sampleVec[0]) - int(desc[i]));
        rowB[0] = rowA[0] + flawCost + computeError(sampleVec[0], thMinVec[i], thMaxVec[i]);
        int cost = min(rowA[1] + flawCost, rowA[0]);
//        rowB[1] = cost + abs(int(sampleVec[1]) - int(desc[i]));
        rowB[1] = cost + computeError(sampleVec[1], thMinVec[i], thMaxVec[i]);
        for (int j = 2; j < sampleVec.size(); j++)
        {
            cost = min(min(rowA[j] + flawCost, rowA[j - 1]), rowA[j - 2] + flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        swap(rowA, rowB);
    }
    vector<int> rowC(sampleVec.size()); //center cost
    swap(rowA, rowC);
    
    //match the second half (from the last pixel to first)
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(sampleVec[i] - desc.back());
        rowA[i] = computeError(sampleVec[i], thMinVec.back(), thMaxVec.back());
    }
    for (int i = desc.size() - 2; i > HALF_LENGTH; i--)
    {
        for (int j = 0; j < sampleVec.size() - 2; j++)
        {
            int cost = min(min(rowA[j] + flawCost, rowA[j + 1]), rowA[j + 2] + flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        int j = sampleVec.size() - 2;
        int cost = min(rowA[j] + flawCost, rowA[j + 1]);
//        rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        
//        rowB.back() = rowA.back() + flawCost + abs(int(sampleVec.back()) - int(desc[i]));
        rowB.back() = rowA.back() + flawCost + computeError(sampleVec.back(), thMinVec[i], thMaxVec[i]);
        swap(rowA, rowB);
    }
    
    //accumulate the cost
    for (int i = 0; i < sampleVec.size() - 2; i++)
    {
        rowC[i] += min(min(rowA[i] + flawCost, rowA[i + 1]), rowA[i + 2] + flawCost);
    }
    int i = rowC.size() - 2;
    rowC[i] += min(rowA[i] + flawCost, rowA[i + 1]);
    rowC.back() += rowA.back() + flawCost;
    return rowC;
} 
