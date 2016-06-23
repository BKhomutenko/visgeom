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
#include <vector>
#include <array>

// Eigen
#include <Eigen/Eigen>


// Eigen data structures
template<typename T>
using Vector2 = Eigen::Matrix<T, 2, 1>;

template<typename T>
using Vector3 = Eigen::Matrix<T, 3, 1>;

template<typename T>
using Matrix3 = Eigen::Matrix<T, 3, 3>;

using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;

// STL containers
using std::vector;
using std::array;

template<typename T>
using Vector3Vec = vector<Vector3<T>>;

using Vector2dVec = vector<Vector2d>;
using Vector3dVec = vector<Vector3d>;
