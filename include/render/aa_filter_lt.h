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

#pragma once

#include "std.h"

/*
A singleton class to store the kernel values for different number of uniform samples
*/
class LookupFilter
{
public:
    static const LookupFilter & instance();
    
    const vector<double> & getKernel(int length) const;
    
    double computeOrigin(int length) const;
    
    double computeStep(int length) const;
    
    double getRadius() const;
    
protected:
    LookupFilter();
    
    vector<vector<double>> _kernelVec;
    
    double computeKernel(double x) const;
    
    const int MAX_SIZE = 25;
    const double SIGMA = 0.55;
    /*
    sigma_x * sigma_f = 1 / (2 * pi) for gaussian filters
    if we want to cut everything that is higher than f = 0.5, then sigma_f must be at least 0.3
    from where we get sigma_x ~= 0.5
    */
};

