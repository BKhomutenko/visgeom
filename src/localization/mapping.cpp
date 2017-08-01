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
SLAM system
*/

#include "localization/mapping.h"


//Speed estimation and extrapolation are to be added
//Assumed that the data arrives in the chronological order

    
void PhotometricMapping::feedOdometry(const Transf & xi)
{    
    //integrate
    _xiLocal = _xiLocal.composeInverse(_xiOdom).compose(xi);
    
    //refresh
    _xiOdom = xi;
}

void PhotometricMapping::feedImage(const Mat8u & img)
{
    _xiOdomImage   
}

void PhotometricMapping::reInit(const Transf & xi)
{
    _state = ST_SELECT;
    _xi = xi;
}

int PhotometricMapping::select()
{
    int res = -1;
    double bestDist = DOUBLE_MAX;
    for (int i = 0; i < _frameVec.size(); i++)
    {
        Transf delta = xi.inverseCompose(_frameVec[i].xi);
        double r = delta.rot().squaredNorm();
        double d = delta.trans().squaredNorm();
        if (r > _params.angualrThreshSq or
            d > _params.distThreshSq) continue;
        if (r + d < bestDist)
        {
            bestDist = r + d;
            res = i;
        }
    }
    return res;
}

void PhotometricMapping::newFrame(const Mat8u & img, int frameIdx)
{
    EnhancedSgm sgm(getCameraMotion().inverse(), _camera, _camera, _sgmParams);
    DepthMap depthNew;
    Mat8u & baseImg = (frameIdx == -1) ? 
                        _tmpFrame.img :
                        _frameVec[frameIdx].img;
                        
    sgm.computeStereo(img, baseImg, depthNew);
}

Transf PhotometricMapping::localizeMI(const Mat8u & img)
{
    
}

Transf PhotometricMapping::localizePhoto(const Mat8u & img)
{

}








