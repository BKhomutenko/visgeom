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

#include "render/aa_filter_lt.h"

const LookupFilter & LookupFilter::instance()
{
    static LookupFilter filter;
    return filter;
}

double LookupFilter::computeKernel(double x) const
{
    return exp(-x*x / (2 * SIGMA * SIGMA)) - 0. * exp(-x*x / (1.5 * SIGMA * SIGMA));
}

const vector<double> & LookupFilter::getKernel(int length) const
{
    assert(length < MAX_SIZE);
    return _kernelVec[length];
}

double LookupFilter::computeOrigin(int length) const
{
    return 3 * SIGMA * (-1 + 1. / length);
}
    
double LookupFilter::computeStep(int length) const
{
    return 6 * SIGMA / length;
}

double LookupFilter::getRadius() const
{
    return 5 * SIGMA;
}

LookupFilter::LookupFilter()
{
    //precompute the integral array
    const int SAMPLE_COUNT = 600;
    
    vector<double> kernelData(SAMPLE_COUNT + 1);
    vector<double> integralKernel(SAMPLE_COUNT + 1);
    integralKernel[0] = 0;
    const double step = computeStep(SAMPLE_COUNT);
    double x = computeOrigin(SAMPLE_COUNT);
    for (int i = 0; i < SAMPLE_COUNT; i++, x += step)
    {
        kernelData[i] = computeKernel(x);
        integralKernel[i + 1] = integralKernel[i] + kernelData[i]; 
    }
    kernelData.back() = 0;
    // normalize the kernel
    double normalizer = 1. / integralKernel.back();
    integralKernel.back() = 1;
    for (int i = 0; i < SAMPLE_COUNT; i++, x += step)
    {
        integralKernel[i] *= normalizer; 
        kernelData[i] *= normalizer; 
    }
    
    _kernelVec.emplace_back(); // 0-size
    
    _kernelVec.emplace_back(); // 1-size
    _kernelVec.back().push_back(1);
    
    _kernelVec.emplace_back(); // 2-size
    _kernelVec.back().push_back(0.5);
    _kernelVec.back().push_back(0.5);
    
    //the rest is generated automatically
    for (int length = 3; length < MAX_SIZE; length++)
    {
        _kernelVec.emplace_back(); // 2-size
        vector<double> & subsampleKernel = _kernelVec.back();
        const double fracSamples = double(SAMPLE_COUNT) / length;
        double prevVal = 0;
        for (int i = 0; i < length; i++)
        {
            double fracIdx = fracSamples * (i + 1);
            int idx = int(fracIdx);
            double surplus = fracIdx - idx;
            double curVal = integralKernel[idx] + surplus * kernelData[idx];
            subsampleKernel.push_back(curVal - prevVal);
            prevVal = curVal;            
        }
    }
}
 
 
    
