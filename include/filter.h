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
A timer to profile functions
*/

#pragma once

#include "std.h"

//TODO bring it away
const std::array<int32_t, 3> KERNEL_3 = {2, 3, 2};
const std::array<int32_t, 5> KERNEL_5 = {2, 4, 5, 4, 2};
const std::array<int32_t, 7> KERNEL_7 = {2, 3, 4, 5, 4, 3, 2};
const std::array<int32_t, 9> KERNEL_9 = {2, 3, 4, 4, 5, 4, 4, 3, 2};


//TODO change filters
const std::array<int32_t, 3> WAVE_3 = {1, -2, 1};
const std::array<int32_t, 5> WAVE_5 = {1, -4, 6, -4, 1};
const std::array<int32_t, 7> WAVE_7 = {1, -6, 15, -20, 15, -6, 1};
const std::array<int32_t, 9> WAVE_9 = {1, -8, 28, -56, 70, -56, 28, -8, 1};

const int NORMALIZER_3 = 7;
const int NORMALIZER_5 = 17;
const int NORMALIZER_7 = 23;
const int NORMALIZER_9 = 31;

const int WAVE_NORM_3 = 2;
const int WAVE_NORM_5 = 8;
const int WAVE_NORM_7 = 30;
const int WAVE_NORM_9 = 90;

template<typename Titer1, typename Titer2, typename Tdata>
Tdata filter(Titer1 kernelIter, Titer1 kernelEnd, Titer2 dataIter, Tdata val)
{
    for (; kernelIter != kernelEnd; ++kernelIter, ++dataIter)
    {
        val += (*kernelIter) * (*dataIter);
    }
    return val;
}
