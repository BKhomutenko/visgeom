
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
Type definitions
*/

#pragma once

// STL
#include "std.h"

// Eigen
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
// Eigen data structures
template<typename T>
using Vector2 = Eigen::Matrix<T, 2, 1>;

template<typename T>
using Vector3 = Eigen::Matrix<T, 3, 1>;

template<typename T>
using Matrix3 = Eigen::Matrix<T, 3, 3>;

using Covector2d = Eigen::Matrix<double, 1, 2>;
using Covector3d = Eigen::Matrix<double, 1, 3>;
using Covector6d = Eigen::Matrix<double, 1, 6>;
using Matrix23drm = Eigen::Matrix<double, 2, 3, Eigen::RowMajor>;
using Matrix32d = Eigen::Matrix<double, 3, 2>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Matrix6drm = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>;
using Matrix56d = Eigen::Matrix<double, 5, 6>;
using Matrix5d = Eigen::Matrix<double, 5, 5>;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Vector5d = Eigen::Matrix<double, 5, 1>;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::Map;

using Eigen::JacobiSVD;

template<typename T>
using Vector3Vec = std::vector<Vector3<T>>;
template<typename T>
using Vector2Vec = std::vector<Vector2<T>>;

using Vector2dVec = std::vector<Vector2d>;
using Vector3dVec = std::vector<Vector3d>;

using Vector2iVec = std::vector<Vector2i>;
using Vector3iVec = std::vector<Vector3i>;

inline Vector2i round(const Vector2d & x)
{
    return Vector2i(round(x[0]), round(x[1]));
}



