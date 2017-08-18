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
STL types and methods
*/

#pragma once

#include <assert.h>
#include <vector>
#include <array>
#include <map>
#include <queue>
#include <list>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>
#include <random>

// Containers
using std::vector;
using std::queue;
using std::list;
using std::array;
using std::pair;
using std::map;
using Array6d = array<double, 6>;
using std::make_pair;

// Algorithms
using std::copy;
using std::sort;
using std::fill;
using std::accumulate;
using std::max_element;
using std::min_element;
using std::distance;

// Math
using std::min;
using std::max;
using std::sin;
using std::cos;
using std::exp;
using std::log;
using std::abs;

// Random
using std::mt19937;

//constants
const double HALF_PI = M_PI / 2;
const double DOUBLE_INF = std::numeric_limits<double>::infinity();
const double DOUBLE_MAX = std::numeric_limits<double>::max();
const double DOUBLE_BIG = 1e15;
const double DOUBLE_SMALL = 1e-6;

inline int sign(double x)
{
    return 2*int(x > 0) - 1;
}


