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
geometry function definition
*/

#pragma once

template<typename T>
inline T sinc(const T & x)
{
    if (x == T(0.))
        return T(1.);
    return sin(x)/x;
}

template<typename T>
Matrix3<T> rotationMatrix(const Vector3<T> & v)
{
    Matrix3<T> R;
    T th = v.norm();
    if (th < 1e-5)
    {
        //Rotational part in case when th is small
        
	    R <<        T(1.),   -v(2),    v(1),
	          v(2),          T(1.),   -v(0),
	         -v(1),    v(0),          T(1.);
    }
    else
    {
        //Rodrigues formula
        T u1 = v(0) / th;
        T u2 = v(1) / th;
        T u3 = v(2) / th;
        T sinth = sin(th);
        T costhVar = T(1.) - cos(th);
        
        R(0, 0) = T(1.) + costhVar * (u1 * u1 - T(1.));
        R(1, 1) = T(1.) + costhVar * (u2 * u2 - T(1.));
        R(2, 2) = T(1.) + costhVar * (u3 * u3 - T(1.));

        R(0, 1) = -sinth*u3 + costhVar * u1 * u2;
        R(0, 2) = sinth*u2 + costhVar * u1 * u3;
        R(1, 2) = -sinth*u1 + costhVar * u2 * u3;

        R(1, 0) = sinth*u3 + costhVar * u2 * u1;
        R(2, 0) = -sinth*u2 + costhVar * u3 * u1;
        R(2, 1) = sinth*u1 + costhVar * u3 * u2;
    }
    return R;
}

template<typename T>
Vector3<T> rotationVector(const Matrix3<T> & R)
{
    return Quaternion<T>(R).toRotationVector;
}

template<typename T>
inline Matrix3<T> hat(const Vector3<T> & u)
{
    Matrix3<T> M;
    M << 0, -u(2), u(1),   u(2), 0, -u(0),   -u(1), u(0), 0;
    return M;
}

template<typename T>
Matrix3<T> interRotOmega(const Vector3<T> & v)
{
    Matrix3<T> B;
    T theta = v.norm();
    if (theta < 1e-5)
    {
        Vector3<T> vHalf = v / 2;
	    B <<        T(1.),    vHalf(2),    -vHalf(1),
	            -vHalf(2),       T(1.),     vHalf(0),
	             vHalf(1),   -vHalf(0),        T(1.);
    }
    else
    {
        Matrix3<T> uhat = hat<T>(v / theta);
        T thetaHalf = theta/2;
        T S1 = sinc(theta);
        T S2 = sinc(thetaHalf);
        T K2 = T(1.) - S1 / (S2 * S2);
        B = Matrix3<T>::Identity() - thetaHalf*uhat + K2*uhat*uhat;
    }
    return B;
}

template<typename T>
Matrix3<T> interOmegaRot(const Vector3<T> & v)
{
    Matrix3<T> B;
    T theta = v.norm();
    if (theta < 1e-5)
    {
        Vector3<T> vHalf = v / T(2.);
	    B <<        T(1.),   -vHalf(2),     vHalf(1),
	             vHalf(2),       T(1.),    -vHalf(0),
	             -vHalf(1),   vHalf(0),        T(1.);
    }
    else
    {
        Matrix3<T> uhat = hat<T>(v / theta);
        T thetaHalf = theta / T(2.);
        T K1 = sinc(thetaHalf);
        K1 = thetaHalf * K1 * K1;
        T K2 = (T(1.) - sinc(theta));
        B = Matrix3<T>::Identity() + K1*uhat + K2*uhat*uhat;
    }
    return B;
}


