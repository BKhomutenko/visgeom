
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
using Eigen::Matrix;

template<typename T>
using Vector2 = Matrix<T, 2, 1>;

template<typename T>
using Vector3 = Matrix<T, 3, 1>;

template<typename T>
using Matrix3 = Matrix<T, 3, 3>;

template<int R, int C>
using Matrixd = Matrix<double, R, C>;

using Eigen::Dynamic;
using Eigen::RowMajor;
using Covector2d = Matrix<double, 1, 2>;
using Covector3d = Matrix<double, 1, 3>;
using Covector6d = Matrix<double, 1, 6>;
using Matrix23drm = Matrix<double, 2, 3, RowMajor>;
using Matrix6d = Matrix<double, 6, 6>;
using Matrix6drm = Matrix<double, 6, 6, RowMajor>;
using Matrix5d = Matrix<double, 5, 5>;
using MatrixXdrm = Matrix<double, Dynamic, Dynamic, RowMajor>;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Vector6d = Matrix<double, 6, 1>;
using Vector5d = Matrix<double, 5, 1>;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::Map;

using Eigen::JacobiSVD;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;

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

template<typename JacobiSvdType, typename MatrixType>
void pinv(JacobiSvdType & svd, MatrixType & pinvmat)
{
    const double PINV_TOLERANCE = 1.e-6; // choose your tolerance wisely!
    auto singularValues = svd.singularValues();
    for (int i = 0; i < singularValues.rows(); i++) {
        if ( singularValues(i) > PINV_TOLERANCE)
        {
            singularValues(i) = 1.0 / singularValues(i);
        }
        else singularValues(i) = 0;
    }
    pinvmat= (svd.matrixV()*singularValues.asDiagonal()*svd.matrixU().transpose());
}

