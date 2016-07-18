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

#include <chrono>

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    
    void reset() 
    { 
        beg_ = clock_::now(); 
    }
    
    double elapsed() const 
    { 
        return std::chrono::duration_cast<second_>(clock_::now() - beg_).count(); 
    }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

