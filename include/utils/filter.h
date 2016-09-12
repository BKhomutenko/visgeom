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
const std::array<int32_t, 3> WAVE_3 = {2, -4, 2};
const std::array<int32_t, 5> WAVE_5 = {1, -3, 4, -3, 1};
const std::array<int32_t, 7> WAVE_7 = {1, -4, 8, -10, 8, -4, 1};
const std::array<int32_t, 9> WAVE_9 = {3, -15, 40, -70, 84, -70, 40, -15, 3};

const int NORMALIZER_3 = 7;
const int NORMALIZER_5 = 17;
const int NORMALIZER_7 = 23;
const int NORMALIZER_9 = 31;

const int WAVE_NORM_3 = 5;
const int WAVE_NORM_5 = 6;
const int WAVE_NORM_7 = 16;
const int WAVE_NORM_9 = 143;

inline int initKernel(vector<int32_t> & outVec, int length)
{
    outVec.resize(length);
    switch (length)
    {
    case 3:
        copy(KERNEL_3.begin(), KERNEL_3.end(), outVec.begin());
        return NORMALIZER_3;
    case 5:
        copy(KERNEL_5.begin(), KERNEL_5.end(), outVec.begin());
        return NORMALIZER_5;
    case 7:
        copy(KERNEL_7.begin(), KERNEL_7.end(), outVec.begin());
        return NORMALIZER_7;
    case 9:
        copy(KERNEL_9.begin(), KERNEL_9.end(), outVec.begin());
        return NORMALIZER_9;
    default:
        throw std::logic_error("length is not in {3, 5, 7, 9}");
    }
}

inline int initWave(vector<int32_t> & outVec, int length)
{
    outVec.resize(length);
    switch (length)
    {
    case 3:
        copy(WAVE_3.begin(), WAVE_3.end(), outVec.begin());
        return WAVE_NORM_3;
    case 5:
        copy(WAVE_5.begin(), WAVE_5.end(), outVec.begin());
        return WAVE_NORM_5;
    case 7:
        copy(WAVE_7.begin(), WAVE_7.end(), outVec.begin());
        return WAVE_NORM_7;
    case 9:
        copy(WAVE_9.begin(), WAVE_9.end(), outVec.begin());
        return WAVE_NORM_9;
    default:
        throw std::logic_error("length is not in {3, 5, 7, 9}");
    }
}

template<typename Titer1, typename Titer2, typename Tdata>
Tdata filter(Titer1 kernelIter, Titer1 kernelEnd, Titer2 dataIter, Tdata val)
{
    for (; kernelIter != kernelEnd; ++kernelIter, ++dataIter)
    {
        val += (*kernelIter) * (*dataIter);
    }
    return val;
}

// bias = (sum(data2) - sum(data1)) / LENGTH
template<typename Titer1, typename Titer2, typename Tdata>
Tdata biasedAbsDiff(Titer1 kernelIter, Titer1 kernelEnd, Titer2 data1Iter,
            Titer2 data2Iter, Tdata bias)
{
    Tdata res = Tdata(0);
    for (; kernelIter != kernelEnd; ++kernelIter, ++data1Iter, ++data2Iter)
    {
        res += (*kernelIter) * abs(*data1Iter - *data2Iter + bias);
    }
    return res;
}

template<typename Titer1, typename Titer2, typename Tdata>
Tdata biasedAbsDiff(Titer1 kernelIter, Titer1 kernelEnd, Titer2 data1Iter,
            Titer2 data2Iter, Tdata bias, int step2)
{
    Tdata res = Tdata(0);
    for (; kernelIter != kernelEnd; ++kernelIter, ++data1Iter, advance(data2Iter, step2))
    {
        res += (*kernelIter) * abs(*data1Iter - *data2Iter + bias);
    }
    return res;
}

template<typename Titer1, typename Tdata>
Tdata totalVariation(Titer1 dataIter, Titer1 dataEnd, Tdata bias)
{
    Tdata res = Tdata(0);
    Tdata sample1 = *dataIter;
    dataIter++;
    for (; dataIter != dataEnd; ++dataIter)
    {
        Tdata sample2 = *dataIter;
        res += abs(sample1 - sample2);
        sample1  = sample2;
    }
    return res;
}

