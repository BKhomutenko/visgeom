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
#include "camera/generic_camera.h"

class StereoEpipoles
{
public:
    StereoEpipoles(const ICamera * camera1, const ICamera * camera2, const Transformation<double> & Transform12)
    {
        if (camera1->projectPoint(Transform12.trans(), epipole1))
        {
            epipoleInverted1 = false;
        }
        else
        {
            camera1->projectPoint(-Transform12.trans(), epipole1);
            epipoleInverted1 = true;
        }
        
        if (camera2->projectPoint(Transform12.transInv(), epipole2))
        {
            epipoleInverted2 = false;
        }
        else
        {
            camera2->projectPoint(-Transform12.transInv(), epipole2);
            epipoleInverted2 = true;
        }
        cout << "Epipoles' inversion " << epipoleInverted1 << " " << epipoleInverted2 << endl;
        
        epipolePx1 = round(epipole1);
        epipolePx2 = round(epipole2);
        cout << "Epipoles : " << epipolePx1.transpose() << " " << epipolePx2.transpose() << endl;
    
    }
    
    const Vector2d & getFirst() const { return epipole1; }
    const Vector2i & getFirstPx() const { return epipolePx1; }
    const Vector2d & getSecond() const { return epipole2; }
    const Vector2i & getSecondPx() const { return epipolePx2; }
    bool firstIsInverted() const { return epipoleInverted1; } 
    bool secondIsInverted() const { return epipoleInverted2; }
    
private:
    bool epipoleInverted1, epipoleInverted2;
    Vector2d epipole1, epipole2;  // projection of the first camera center onto the second camera
    Vector2i epipolePx1, epipolePx2;
};

