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
All necessary geometric transformations.
This header-based package is completely independent of the rest of the project
and can be used as a stand-alone lignt-weight 3D geometry library.
Requires Eigen3
*/

#pragma once

//Eigen
#include <Eigen/Eigen>

#define ZERO T(0.)

#define T_PI T(M_PI)

template<typename T>
using Vector2 = Eigen::Matrix<T, 2, 1>;
template<typename T>
using Vector3 = Eigen::Matrix<T, 3, 1>;
template<typename T>
using Matrix3 = Eigen::Matrix<T, 3, 3>;
template<typename T>
using Vector3Vec = std::vector<Vector3<T>>;

template<typename T>
T sinc(const T & x);

template<typename T>
T normalizeAngle(const T & th);

template<typename T>
Matrix3<T> rotationMatrix(const Vector3<T> & v);

template<typename T>
Vector3<T> rotationVector(const Matrix3<T> & R);

template<typename T>
Matrix3<T> hat(const Vector3<T> & u);

template<typename T>
Matrix3<T> interRotOmega(const Vector3<T> & v);

template<typename T>
Matrix3<T> interOmegaRot(const Vector3<T> & v);

#include "geometry/quaternion.h"
#include "geometry/transformation.h"
#include "geometry/geometry_core.h"


