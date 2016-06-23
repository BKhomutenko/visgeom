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
Template-based quaternion implementation
*/

#pragma once

template<typename T>
class Quaternion
{
public:
    Quaternion() {}
    Quaternion(T x, T y, T z, T w) : data{x, y, z, w} {}
    Quaternion(const Vector3<T> & rot)
    {
        T theta = rot.norm();
        if ( abs(theta) < 1e-6 )
        {
            x = rot(0) / 2.;
            y = rot(1) / 2.;
            z = rot(2) / 2.;
            w = T(1.);  
        }
        else
        {
            Vector3<T> u = rot / theta;
            T s = sin(theta / 2.);
            x = u(0) * s;
            y = u(1) * s;
            z = u(2) * s;
            w = cos(theta / 2.); 
        }
    } 
    
    Quaternion(const Matrix3<T> & R)
    {
	    w = sqrt(1.0 + R.trace()) / 2.0;
	    T w4 = (4.0 * w);
	    x = (R(2, 1) - R(1, 2)) / w4 ;
	    y = (R(0, 2) - R(2, 0)) / w4 ;
	    z = (R(1, 0) - R(0, 1)) / w4 ;
    }
    
    Vector3<T> rotate(const Vector3<T> & v) const
    {
        T t1 =   w*x;
        T t2 =   w*y;
        T t3 =   w*z;
        T t4 =  -x*x;
        T t5 =   x*y;
        T t6 =   x*z;
        T t7 =  -y*y;
        T t8 =   y*z;
        T t9 =  -z*z;
        
        const T & v1 = v(0);
        const T & v2 = v(1);
        const T & v3 = v(2);
        
        T v1new = 2.*((t7 + t9)*v1 + (t5 - t3)*v2 + (t2 + t6)*v3 ) + v1;
        T v2new = 2.*((t3 + t5)*v1 + (t4 + t9)*v2 + (t8 - t1)*v3 ) + v2;
        T v3new = 2.*((t6 - t2)*v1 + (t1 + t8)*v2 + (t4 + t7)*v3 ) + v3;
        
        return Vector3<T>(v1new, v2new, v3new);
    } 
    
    Vector3<T> toRotationVector() const
    {
        T s = sqrt(x*x + y*y + z*z);
        Vector3<T> u(x, y, z);
        if (s < 1e-5)
        {
            return u * T(2.);
        }
        else
        {
            T th = 2. * atan2(s, w);
            return u / s * th;
        }
    } 
    
    Quaternion inv() const
    {
        return Quaternion(-x, -y, -z, w);
    }
    
    Quaternion operator*(const Quaternion & q) const
    {
        const T & x2 = q.x;
        const T & y2 = q.y;
        const T & z2 = q.z;
        const T & w2 = q.w;
        
        T wn = w*w2 - x*x2 - y*y2 - z*z2;
        T xn = w*x2 + x*w2 + y*z2 - z*y2;
        T yn = w*y2 - x*z2 + y*w2 + z*x2;
        T zn = w*z2 + x*y2 - y*x2 + z*w2;
        
        return Quaternion(xn, yn, zn, wn);
    } 

    friend std::ostream& operator << (std::ostream & os, const Quaternion & Q)
    {
        os << Q.x << " " << Q.y << " " << Q.z << " " << Q.w;
        return os;
    }
private:
    T data[4];
    T & x = data[0];
    T & y = data[1];
    T & z = data[2];
    T & w = data[3];
};

