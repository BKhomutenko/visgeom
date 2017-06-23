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
#include "reconstruction/stereo_misc.h"

enum EpipoleResult : uint32_t {EPIPOLE_INVERTED = 1, EPIPOLE_TOO_CLOSE = 2};

class StereoEpipoles
{
public:

    StereoEpipoles() {}
    
    StereoEpipoles(const ICamera * camera1, const ICamera * camera2,
            const Transf & Transform12);
    
    const Vector2d & get(CameraIdx idx, uint32_t result) const 
    { 
        if (result & EPIPOLE_INVERTED) return antiEpipole[idx];
        else return epipole[idx];
    }
    
    const Vector2d & get(CameraIdx idx) const 
    { 
        if (epipoleProjected[idx]) return epipole[idx];
        else return antiEpipole[idx];
    }
    
    const Vector2i & getPx(CameraIdx idx, uint32_t result) const 
    { 
        if (result & EPIPOLE_INVERTED) return antiEpipolePx[idx];
        else return epipolePx[idx];
    }
    
    const Vector2i & getPx(CameraIdx idx) const 
    { 
        if (epipoleProjected[idx]) return epipolePx[idx];
        else return antiEpipolePx[idx];
    }
    
//    const Vector2d & getSecond(EpipoleResult result) const 
//    { 
//        if (result & EPIPOLE_INVERTED) return antiEpipole2;
//        else return epipole2;
//    }
//    
//    const Vector2d & getSecond() const 
//    { 
//        if (epipole2projected) return epipole2;
//        else return antiEpipole2;
//    }
//    
//    const Vector2i & getSecondPx(EpipoleResult result) const 
//    { 
//        if (result & EPIPOLE_INVERTED) return antiEpipolePx2;
//        else return epipolePx2;
//    }
//    
//    const Vector2i & getSecondPx() const 
//    { 
//        if (epipole2projected) return epipolePx2;
//        else return antiEpipolePx2;
//    }
    
    //TODO testing
    uint32_t chooseEpipole(CameraIdx idx, const Vector2i pt, int threshSquared = 0) const;
    
private:
    array<bool, 2> epipoleProjected;
    array<bool, 2> antiEpipoleProjected;
    array<Vector2d, 2> epipole;  // projection of the first camera center onto the second camera
    array<Vector2d, 2> antiEpipole;  // projection of the first camera center onto the second camera
    array<Vector2i, 2> epipolePx;
    array<Vector2i, 2> antiEpipolePx;
};

