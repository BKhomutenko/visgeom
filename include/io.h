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
input-outpuit
*/

#pragma once

// STL
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using std::string;

using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::endl;
using std::setw;
using std::flush;
using std::to_string;


template <typename T>
void printVector(const T begin, const T end, int w = 8)
{
    for (auto it = begin; it != end; it++)
    {
        cout << setw(w) << *it;
    }
    cout << endl;
}

template <typename T>
void printPointVector(const T begin, const T end)
{
    for (auto it = begin; it != end; it++)
    {
        cout << (*it).transpose() << endl;
    }
}

