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


#include "reconstruction/epipoles.h"

    
StereoEpipoles::StereoEpipoles(const ICamera * camera1, const ICamera * camera2,
        const Transf & Transform12)
{
    epipole1projected = camera1->projectPoint(Transform12.trans(), epipole1);
    antiEpipole1projected = camera1->projectPoint(-Transform12.trans(), antiEpipole1); 
    epipole2projected = camera2->projectPoint(Transform12.transInv(), epipole2);
    antiEpipole2projected =  camera2->projectPoint(-Transform12.transInv(), antiEpipole2);
    
    
    //bound the epipole coordinates
    vector<Vector2d *> epipoleVec = {&epipole1, &epipole2, &antiEpipole1, &antiEpipole2};
    for (auto xptr : epipoleVec)
    {
        Vector2d & x = *xptr;
        for (int i = 0; i < 2; i++)
        {
            if (x[i] > 1e6) x[i] = 1e6;
            if (x[i] < -1e6) x[i] = -1e6; 
        }
        
    }
    
    
    if (epipole1projected) epipolePx1 = round(epipole1);
    if (antiEpipole1projected) antiEpipolePx1 = round(antiEpipole1);
    if (epipole2projected) epipolePx2 = round(epipole2);
    if (antiEpipole2projected) antiEpipolePx2 = round(antiEpipole2);
}
    

    
//TODO testing
bool StereoEpipoles::useInvertedEpipoleSecond(const Vector2i pt) const
{
    if (epipole2projected and antiEpipole2projected)
    {
        double dist = (pt - epipolePx2).squaredNorm();
        double antiDist = (pt - antiEpipolePx2).squaredNorm();
        return dist > antiDist;
    }
    else if (epipole2projected) return false;
    else if (antiEpipole2projected) return true;
    else
    {
        throw;
    }
    
}

bool StereoEpipoles::useInvertedEpipoleFirst(const Vector2i pt) const
{
    if (epipole1projected and antiEpipole1projected)
    {
        double dist = (pt - epipolePx1).squaredNorm();
        double antiDist = (pt - antiEpipolePx1).squaredNorm();
        return dist > antiDist;
    }
    else if (epipole1projected) return false;
    else if (antiEpipole1projected) return true;
    else
    {
        throw;
    }
    
}

