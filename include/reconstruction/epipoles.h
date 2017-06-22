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
A utility class for stereo reconstruction
Computes and keeps the epipoles or aniepipoles
*/

#pragma once

#include "io.h"

#include "eigen.h"
#include "geometry/geometry.h"
#include "projection/generic_camera.h"

class StereoEpipoles
{
public:

    StereoEpipoles() {}
    
    StereoEpipoles(const ICamera * camera1, const ICamera * camera2,
            const Transf & Transform12);
    
    const Vector2d & getFirst(bool inverted) const 
    { 
        if (inverted) return antiEpipole1;
        else return epipole1;
    }
    
    const Vector2d & getFirst() const 
    { 
        if (epipole1projected) return epipole1;
        else return antiEpipole1;
    }
    
    const Vector2i & getFirstPx(bool inverted) const 
    { 
        if (inverted) return antiEpipolePx1;
        else return epipolePx1;
    }
    
    const Vector2i & getFirstPx() const 
    { 
        if (epipole1projected) return epipolePx1;
        else return antiEpipolePx1;
    }
    
    const Vector2d & getSecond(bool inverted) const 
    { 
        if (inverted) return antiEpipole2;
        else return epipole2;
    }
    
    const Vector2d & getSecond() const 
    { 
        if (epipole2projected) return epipole2;
        else return antiEpipole2;
    }
    
    const Vector2i & getSecondPx(bool inverted) const 
    { 
        if (inverted) return antiEpipolePx2;
        else return epipolePx2;
    }
    
    const Vector2i & getSecondPx() const 
    { 
        if (epipole2projected) return epipolePx2;
        else return antiEpipolePx2;
    }
    
    //TODO testing
    bool useInvertedEpipoleSecond(const Vector2i pt) const;
    bool useInvertedEpipoleFirst(const Vector2i pt) const;
    
//private:
    bool epipole1projected, epipole2projected;
    bool antiEpipole1projected, antiEpipole2projected;
    Vector2d epipole1, epipole2;  // projection of the first camera center onto the second camera
    Vector2d antiEpipole1, antiEpipole2;  // projection of the first camera center onto the second camera
    Vector2i epipolePx1, epipolePx2;
    Vector2i antiEpipolePx1, antiEpipolePx2;
};

