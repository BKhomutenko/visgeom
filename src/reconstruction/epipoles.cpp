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

void limitVector(Vector2d & x)
{
    const double MAX_VAL = 1e6;
    if (x[0] > MAX_VAL) x[0] = MAX_VAL;
    if (x[0] < -MAX_VAL) x[0] = -MAX_VAL;
    if (x[1] > MAX_VAL) x[1] = MAX_VAL;
    if (x[1] < -MAX_VAL) x[1] = -MAX_VAL;
}
    
StereoEpipoles::StereoEpipoles(const ICamera * camera1, const ICamera * camera2,
        const Transf & Transform12)
{
    epipoleProjected[CAMERA_1] = camera1->projectPoint(Transform12.trans(), epipole[CAMERA_1]);
    antiEpipoleProjected[CAMERA_1] = camera1->projectPoint(-Transform12.trans(), antiEpipole[CAMERA_1]); 
    
    Vector3d transInv = Transform12.transInv();
    epipoleProjected[CAMERA_2] = camera2->projectPoint(transInv, epipole[CAMERA_2]);
    antiEpipoleProjected[CAMERA_2] =  camera2->projectPoint(-transInv, antiEpipole[CAMERA_2]);
    
    
    
    for (int i = 0; i < 2; i++)
    {
        //bound the epipole coordinates
        limitVector(epipole[i]);
        limitVector(antiEpipole[i]);
        
        // convert to integer pixels
        if (epipoleProjected[i]) epipolePx[i] = round(epipole[i]);
        if (antiEpipoleProjected[i]) antiEpipolePx[i] = round(antiEpipole[i]);
    }
}
    

    
//TODO testing
uint32_t StereoEpipoles::chooseEpipole(CameraIdx idx, const Vector2i pt, int threshSquared) const
{
    uint32_t res = 0;
    if (epipoleProjected[idx] and antiEpipoleProjected[idx])
    {
        double dist = (pt - epipolePx[idx]).squaredNorm();
        double antiDist = (pt - antiEpipolePx[idx]).squaredNorm();
        if (antiDist < dist)
        {
            res |= EPIPOLE_INVERTED;
            if (antiDist < threshSquared) res |= EPIPOLE_TOO_CLOSE;
        }
        else if (dist < threshSquared) res |= EPIPOLE_TOO_CLOSE;
        
    }
    else if (epipoleProjected[idx])
    {
        double dist = (pt - epipolePx[idx]).squaredNorm();
        if (dist < threshSquared) res |= EPIPOLE_TOO_CLOSE;
    }
    else if (antiEpipoleProjected[idx])
    {
        double antiDist = (pt - antiEpipolePx[idx]).squaredNorm();
        res |= EPIPOLE_INVERTED;
        if (antiDist < threshSquared) res |= EPIPOLE_TOO_CLOSE;
    }
    else
    {
        throw;
    }
    return res;
}


