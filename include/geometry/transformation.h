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
        os << transfo.mtrans.transpose() << " " << transfo.mrot.transpose();
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
    
    template<typename T2>
    Transformation<T2> cast() const
    {
        Vector3<T2> trans = mtrans.template cast<T2>();
        Vector3<T2> rot = mrot.template cast<T2>();
        return Transformation<T2>(trans, rot);
    }
    
    //TODO check
    Matrix<T, 6, 6> screwTransf() const
    {
        Matrix3<T> R = rotMat();
        Matrix<T, 6, 6> res;
        res.template topLeftCorner<3, 3>() = R;
        res.template topRightCorner<3, 3>() = hat(mtrans) * R;
        res.template bottomLeftCorner<3, 3>() = Matrix3<T>::Zero();
        res.template bottomRightCorner<3, 3>() = R;
        return res;
    }
    
    Matrix<T, 6, 6> screwTransfInv() const
    {
        Matrix3<T> R = rotMatInv();
        Matrix<T, 6, 6> res;
        res.template topLeftCorner<3, 3>() = R;
        res.template topRightCorner<3, 3>() = -R * hat(mtrans);
        res.template bottomLeftCorner<3, 3>() = Matrix3<T>::Zero();
        res.template bottomRightCorner<3, 3>() = R;
        return res;
    }
    
private:
    Vector3<T> mrot;
    Vector3<T> mtrans;

};

using Transf = Transformation<double>;


