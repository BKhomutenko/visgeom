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
#include "io.h"

#include "eigen.h"

// Non-redundant transformation representation
// using translation and angle-axis

template<typename T>
class Transformation
{
public:   
    Transformation() : mrot(ZERO, ZERO, ZERO), mtrans(ZERO, ZERO, ZERO) {}

    Transformation(Vector3<T> trans, Vector3<T> rot) : mtrans(trans), mrot(rot) {}

    Transformation(const Vector3<T> & trans, const Quaternion<T> & qrot)
    : mtrans(trans), mrot(qrot.toRotationVector()) {}
    
    Transformation(const Vector3<T> & trans, const Matrix3<T> & R)
    : mtrans(trans), mrot(rotationVector(R)) {}
    
    Transformation(const T * const data) : mtrans(data), mrot(data + 3) {}

    Transformation(T x, T y, T z, T rx, T ry, T rz)
    : mtrans(x, y, z), mrot(rx, ry, rz) {}

    Transformation(T x, T y, T z, T qx, T qy, T qz, T qw)
    : mtrans(x, y, z), mrot(Quaternion<T>(qx, qy, qz, qw).toRotationVector()) {}

    

    void setParam(Vector3<T> & newTrans, Vector3<T> & newRot)
    {
        mtrans = newTrans;
        mrot = newRot;
    }

    void toRotTrans(Matrix3<T> & R, Vector3<T> & t) const
    {
        t = mtrans;
        R = rotMat();
    }

    void toRotTransInv(Matrix3<T> & R, Vector3<T> & t) const
    {
        R = rotMatInv();
        t = -R*mtrans;
    }

    void scale(const T lambda)
    {
        mrot *= lambda;
        mtrans *= lambda;
    }
    
    Transformation compose(const Transformation & transfo) const
    {
        Transformation res;
        Quaternion<T> q1(mrot), q2(transfo.mrot);
        res.mtrans = q1.rotate(transfo.mtrans) + mtrans;
        Quaternion<T> qres = q1 * q2;
        res.mrot = qres.toRotationVector();
        return res;
    }

    Transformation inverseCompose(const Transformation & transfo) const
    {
        Transformation res;
        Quaternion<T> q1(mrot), q2(transfo.mrot);
        Quaternion<T> q1inv = q1.inv();
        res.mtrans = q1inv.rotate(transfo.mtrans - mtrans);
        Quaternion<T> qres = q1inv * q2;
        res.mrot = qres.toRotationVector();
        return res;
    }

    Transformation composeInverse(const Transformation & transfo) const
    {
        Transformation res;
        Quaternion<T> q1(mrot), q2(transfo.mrot);
        Quaternion<T> q2inv = q2.inv();
        Quaternion<T> qres = q1 * q2inv;
        res.mtrans = mtrans - qres.rotate(transfo.mtrans);
        res.mrot = qres.toRotationVector();
        return res;
    }
    
    Transformation inverse() const
    {
        Transformation res;
        Matrix3<T> R = rotMatInv();
        res.mtrans = -R*mtrans;
        res.mrot = -mrot;
        return res;
    }
    
    const Vector3<T> & trans() const { return mtrans; }

    const Vector3<T> & rot() const { return mrot; }

    Vector3<T> & trans() { return mtrans; }

    Vector3<T> & rot() { return mrot; }

    Quaternion<T> rotQuat() const { return Quaternion<T>(mrot); }

    Matrix3<T> rotMat() const { return rotationMatrix<T>(mrot); }
    Matrix3<T> rotMatInv() const { return rotationMatrix<T>(-mrot); }
    Vector3<T> transInv() const { return -rotMatInv() * mtrans; }
    
    T * rotData() { return mrot.data(); }
    T * transData() { return mtrans.data(); }
    
    const T * rotData() const { return mrot.data(); }
    const T * transData() const { return mtrans.data(); }
    
    friend std::ostream& operator << (std::ostream & os, const Transformation & transfo)
    {
        os << transfo.mtrans.transpose() << " # " << transfo.mrot.transpose();
        return os;
    }

    void transform(const Vector3Vec<T> & src, Vector3Vec<T> & dst) const
    {
        dst.resize(src.size());
        rotate(src, dst);
        for (unsigned int i = 0; i < dst.size(); i++)
        {
            dst[i] = dst[i] + mtrans;
        }
    }

    void inverseTransform(const Vector3Vec<T> & src, Vector3Vec<T> & dst) const
    {
        dst.resize(src.size());
        for (unsigned int i = 0; i < src.size(); i++)
        {
            dst[i] = src[i] - mtrans;
        }
        inverseRotate(dst, dst);
    }

    void transform(const Vector3<T> & src, Vector3<T> & dst) const
    {
        Matrix3<T> R = rotMat();
        dst = R * src + mtrans;
    }
    
    void inverseTransform(const Vector3<T> & src, Vector3<T> & dst) const
    {
        dst = src - mtrans;
        Matrix3<T> R = rotMatInv();
        dst = R * dst;
    }
    
    void rotate(const Vector3Vec<T> & src, Vector3Vec<T> & dst) const
    {
        dst.resize(src.size());
        Matrix3<T> R = rotMat();
        for (unsigned int i = 0; i < src.size(); i++)
        {
            dst[i] = R * src[i];
        }
    }

    void inverseRotate(const Vector3Vec<T> & src, Vector3Vec<T> & dst) const
    {
        dst.resize(src.size());
        Matrix3<T> R = rotMatInv();
        for (unsigned int i = 0; i < src.size(); i++)
        {
            dst[i] = R * src[i];
        }
    }

    std::array<T, 6> toArray() const
    {
        std::array<T, 6> res;
        std::copy(transData(), transData() + 3, res.data());
        std::copy(rotData(), rotData() + 3, res.data() + 3);
        return res;
    }
    
    void toArray(T * const data) const
    {
        std::copy(transData(), transData() + 3, data);
        std::copy(rotData(), rotData() + 3, data + 3);
    }
    
    //TODO put elsewhere
    //jac is a kinematic jacobian  l_dot = jac  *   (  v  )
    //                                              (omega)
    //Memory must be allocated
    void triangulate(const Vector3Vec<T> & xVec1, const Vector3Vec<T> & xVec2,
            T * res1, T * res2 = NULL,
            T * jac1 = NULL, T * jac2 = NULL) //row-major
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
        Matrix3<T> R = rotMat();
        for (int i = 0; i < xVec1.size(); i++)
        {
            const Vector3<T> & p = xVec1[i];
            Vector3<T> q = R * xVec2[i];
            T pq = p.dot(q);
            T pp = p.dot(p);
            T qq = q.dot(q);
            T tp = mtrans.dot(p);
            T tq = mtrans.dot(q);
            T ppqq = pp * qq;
            T delta = -pp * qq + pq * pq;
            if (abs(delta) < ppqq * T(1e-6)) // max 100m for 10cm base
            {
                if (res1 != NULL) res1[i] = ZERO;
                if (res2 != NULL) res2[i] = ZERO;
                if (jac1 != NULL)
                {
                    fill(jac1 + i * 6, jac1 + (i + 1) * 6, ZERO); 
                }
                if (jac2 != NULL)
                {
                    fill(jac2 + i * 6, jac2 + (i + 1) * 6, ZERO); 
                }
            }
            else
            {
                T deltaInv = T(1.) / delta;
                T delta1, delta2;
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
                    Matrix3<T> qSkew = hat(q);
                    Vector3<T> jacDeltaOmega = 2 * pq * qSkew * p;
                    if (jac1 != NULL)
                    {
                        Vector3<T> jacDelta1V = q*pq - p*qq;
                        Vector3<T> jacDelta1Omega = qSkew * (pq * mtrans + tq * p);
                        Map<Vector3<T>> jacLambda1V(jac1 + i * 6);
                        Map<Vector3<T>> jacLambda1Omega(jac1 + i * 6 + 3);
                        jacLambda1V = deltaInv * jacDelta1V;
                        jacLambda1Omega = deltaInv * (jacDelta1Omega - res1[i] * jacDeltaOmega);
                    }
                    if (jac2 != NULL)
                    {
                        Vector3<T> jacDelta2V = q*pp - p*pq;
                        Vector3<T> jacDelta2Omega = qSkew * (pp * mtrans - tp * p);
                        Map<Vector3<T>> jacLambda2V(jac2 + i * 6);
                        Map<Vector3<T>> jacLambda2Omega(jac2 + i * 6 + 3);
                        jacLambda2V = deltaInv * jacDelta2V;
                        jacLambda2Omega = deltaInv * (jacDelta2Omega - res2[i] * jacDeltaOmega);
                    }
                }
            }
        }
    }
    
    // num / denom if denom is big enough
    // linear extrapolation otherwise
    T regDiv(T num, T denom, T EPS)
    {
        if (denom > EPS * num)
        {
            return num / denom;
        }
        else if (num == 0) // singularity, very unlikely
        {
            return T(2.) / EPS;
        }
        else
        {
            return T(2.) / EPS - denom / (num * EPS * EPS);
        }
    }
    
    //below a certain threshold division by delta is replaced by a linear function
    void triangulateRegular(const Vector3Vec<T> & xVec1, const Vector3Vec<T> & xVec2,
            T * res1, T * res2 = NULL,
            T * jac1 = NULL, T * jac2 = NULL) //row-major
    {
        const T EPS(1. / 33.); // max 33m
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
        Matrix3<T> R = rotMat();
        for (int i = 0; i < xVec1.size(); i++)
        {
            const Vector3<T> & p = xVec1[i];
            Vector3<T> q = R * xVec2[i];
            Vector3<T> r = p + q;
            T pq = p.dot(q);
            T pp = p.dot(p);
            T qq = q.dot(q);
            T tp = mtrans.dot(p);
            T tq = mtrans.dot(q);
            T tr = mtrans.dot(r);
            T tt = mtrans.dot(mtrans);
            T rp = r.dot(p);
            T rq = r.dot(q);
            T delta = tp*rq - tq*rp;
            T delta1, delta2;
            if (res1 != NULL) 
            {
                delta1 = tt * rq - tr * tq;
                res1[i] = regDiv(delta1, delta, EPS); // max 33m
//                if (res1[i] != res1[i])
//                { 
//                    std::cout << "OLOLO " << delta << " " << delta1 << std::endl;
//                    cout << *this << endl;
//                    cout << pp << " " << pq << " " << qq << " " << tp << " " << tq << endl;
//                    cout << p.transpose() << "     " << q.transpose() << "   " << mtrans.transpose() << endl;
//                    cout <<  xVec2[i].transpose() << endl << R << endl;
//                }
            }
            if (res2 != NULL)
            {
                delta2 = tt * rp - tr * tp;
                res2[i] = regDiv(delta2, delta, EPS); // max 33m
            }
            
            if (jac1 != NULL or jac2 != NULL)
            {
                T deltaInv = T(1.) / delta;
                Matrix3<T> qSkew = hat(q);
                Vector3<T> jacDeltaV = rq * p - rp * q;
                Vector3<T> jacDeltaOmega = qSkew * (tp * p - tq * p - rp * mtrans);
                if (jac1 != NULL)
                {
                    Vector3<T> jacDelta1V = 2*rq * mtrans - tq * r - tr * q;
                    Vector3<T> jacDelta1Omega = qSkew * (tt * p - 2 * tq * mtrans);
                    Map<Vector3<T>> jacLambda1V(jac1 + i * 6);
                    Map<Vector3<T>> jacLambda1Omega(jac1 + i * 6 + 3);
                    if (delta > EPS * delta1)
                    {
                        jacLambda1V = deltaInv * (jacDelta1V - res1[i] * jacDeltaV);
                        jacLambda1Omega = deltaInv * (jacDelta1Omega - res1[i] * jacDeltaOmega);
                    }
                    else if (delta1 == 0) //very unlikely
                    {
                        jacLambda1V = Vector3<T>(0, 0, 0);
                        jacLambda1Omega = Vector3<T>(0, 0, 0);
                    }
                    else
                    {
                        const T coef = -T(1.) / (EPS * EPS);
                        T deltaInv1 = 1. / delta1;
                        T k = coef * deltaInv1;
                        T k1 = coef * delta * deltaInv1 * deltaInv1;
                        jacLambda1V = k * jacDeltaV - k1 * jacDelta1V;
                        jacLambda1Omega = k * jacDeltaOmega - k1 * jacDelta1Omega;
                    }
                }
                if (jac2 != NULL)
                {
                    Vector3<T> jacDelta2V = 2*rp * mtrans - tp * r - tr * p;
                    Vector3<T> jacDelta2Omega = qSkew * (tt * p - 2 * tp * mtrans);
                    Map<Vector3<T>> jacLambda2V(jac2 + i * 6);
                    Map<Vector3<T>> jacLambda2Omega(jac2 + i * 6 + 3);
                    if (delta > EPS * delta2)
                    {
                        jacLambda2V = deltaInv * (jacDelta2V - res2[i] * jacDeltaV);
                        jacLambda2Omega = deltaInv * (jacDelta2Omega - res2[i] * jacDeltaOmega);
                    }
                    else if (delta2 == 0) //very unlikely
                    {
                        jacLambda2V = Vector3<T>(0, 0, 0);
                        jacLambda2Omega = Vector3<T>(0, 0, 0);
                    }
                    else
                    {
                        const T coef = -T(1.) / (EPS * EPS);
                        T deltaInv2 = 1. / delta2;
                        T k = coef * deltaInv2;
                        T k2 = coef * delta * deltaInv2 * deltaInv2;
                        jacLambda2V = k * jacDeltaV - k2 * jacDelta2V;
                        jacLambda2Omega = k * jacDeltaOmega - k2 * jacDelta2Omega;
                    }
                }
            }
        }
    }
    
    template<typename T2>
    Transformation<T2> cast() const
    {
        Vector3<T2> trans = mtrans.template cast<T2>();
        Vector3<T2> rot = mrot.template cast<T2>();
        return Transformation<T2>(trans, rot);
    }
    
private:
    Vector3<T> mrot;
    Vector3<T> mtrans;

};

using Transf = Transformation<double>;


