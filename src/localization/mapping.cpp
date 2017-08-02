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
PhotometricMapping::PhotometricMapping(const ptree & params):
//    _params(params.get_child("mapping_parameters")),
    _xiBaseCam( readTransform(params.get_child("xi_base_camera")) ),
    _sgmParams(params.get_child("stereo_parameters")),
    _camera( new EnhancedCamera(readVector<double>(params.get_child("camera_params")).data()) ),
    _sparseOdom(_camera, _xiBaseCam),
    _motionStereo(_camera, _camera, params.get_child("stereo_parameters")),
    _odomInit(false),
    _localizer(5, _camera),
    _xiLocal(0, 0, 0, 0, 0, 0),
    _state(MAP_BEGIN)
{
    _localizer.setVerbosity(0);
    _localizer.setXiBaseCam(_xiBaseCam);
}
    
void PhotometricMapping::feedOdometry(const Transf & xi)
{
    if (_odomInit)
    {    //integrate
        _xiLocal = _xiLocal.composeInverse(_xiOdom).compose(xi);
    }
    //refresh
    _xiOdom = xi;
    _odomInit = true;
}

void PhotometricMapping::feedImage(const Mat8u & img)
{
    bool interDistanceOk, mapDistanceOk;
    Transf xiMapFr;
    switch (_state)
    {
    case MAP_BEGIN:
        //store needed data for the initialization;
        _interFrame.xi = _xiLocal;
        img.copyTo(_interFrame.img);
        _sparseOdom.feedData(img, _xiLocal);
        _state = MAP_INIT;
        break;
    case MAP_INIT:
        cout << "INIT" << endl;
        //check whether the init distance is enough;
        if (_xiLocal.trans().norm() < _params.minInitDist) break;
        _sparseOdom.feedData(img, _xiLocal);
        _xiLocal = _sparseOdom.getIncrement();        
        pushInterFrame(img);
        
        _mapIdx = selectMapFrame(_interFrame.xi);
        
        if (_mapIdx == -1)
        {
            _state = MAP_SLAM; //the frame will be put into the map
        }
        else
        {
            _state = MAP_LOCALIZE;
            localizeMI();
        }
        break;
    case MAP_LOCALIZE:
        cout << "LOCALIZE" << endl;
        //The area has been already mapped
        //localize wrt the intermediate frame
        localizePhoto(img);

        //check the constraints
        
        //distance to the inter frame
        interDistanceOk = checkDistance(_xiLocal);
        
        //distane to the map frame
        xiMapFr = _frameVec[_mapIdx].xi.inverseCompose(_interFrame.xi);
        mapDistanceOk = checkDistance(xiMapFr.compose(_xiLocal));
        
        if (not mapDistanceOk) //need to change the map frame
        {
            pushInterFrame(img);
            _mapIdx = selectMapFrame(_interFrame.xi);
            if (_mapIdx != -1) 
            {
                localizeMI();
            }
            else
            {
                _state = MAP_SLAM;
            }
        }
        else if (not interDistanceOk)
        {
            pushInterFrame(img);
            localizeMI();
        }
        else
        {
            _depth = _motionStereo.compute(getCameraMotion(_xiLocal), img, _depth);
            _depth.filterNoise();
        }
        break;
    case MAP_SLAM:
        cout << "SLAM" << endl;
        ///////localize with respect to the intermediate frame/////
        localizePhoto(img);
        
        
        //////check the constraints/////
        if (not checkDistance(_xiLocal))
        {
            pushInterFrame(img);
            _mapIdx = selectMapFrame(_interFrame.xi);
            if (_mapIdx != -1)
            { 
                _state = MAP_LOCALIZE;
                localizeMI();
            }
        }
        else
        {
            _depth = _motionStereo.compute(getCameraMotion(_xiLocal), img, _depth);
            _depth.filterNoise();
        }
        //if the inter frame is too far, push it into the map and 
        //init another intermediate frame
        break;
    }
}

void PhotometricMapping::pushInterFrame(const Mat8u & img)
{
    if (_state == MAP_SLAM)
    {
        _frameVec.emplace_back();
        _interFrame.img.copyTo(_frameVec.back().img);
        _frameVec.back().xi = _interFrame.xi;
    }
    DepthMap newDepth;
    EnhancedSgm sgm(getCameraMotion(_xiLocal).inverse(), _camera, _camera, _sgmParams);
    sgm.computeStereo(img, _interFrame.img, newDepth);
    newDepth.filterNoise();
    if (_state == MAP_INIT)
    {
        _depth = newDepth;
    }
    else
    {
        //project forward the prior depth
        _depth = _depth.wrapDepth(_xiLocal);
        _depth.merge(newDepth);
    }
    img.copyTo(_interFrame.img);
    _interFrame.xi = _interFrame.xi.compose(_xiLocal);
    _xiLocal = Transf(0, 0, 0, 0, 0, 0);
    
    _localizer.setBaseImage(img);
    _motionStereo.setBaseImage(img);
    
}

Transf PhotometricMapping::getCameraMotion(const Transf & xi) const
{
    return _xiBaseCam.inverseCompose(xi).compose(_xiBaseCam);
}

void PhotometricMapping::reInit(const Transf & xi)
{
    if (_state == MAP_SLAM)
    {
        _frameVec.emplace_back();
        _interFrame.img.copyTo(_frameVec.back().img);
        _frameVec.back().xi = _interFrame.xi;
    }
    _state = MAP_BEGIN;
    _xiLocal = xi;
    _odomInit = false;
}

int PhotometricMapping::selectMapFrame(const Transf & xi)
{
    int res = -1;
    double bestDist = DOUBLE_MAX;
    for (int i = 0; i < _frameVec.size(); i++)
    {
        Transf delta = xi.inverseCompose(_frameVec[i].xi);
        double r = delta.rot().squaredNorm();
        double d = delta.trans().squaredNorm();
        if (r > _params.angularThreshSq or
            d > _params.distThreshSq) continue;
        if (r + d < bestDist)
        {
            bestDist = r + d;
            res = i;
        }
    }
    return res;
}

bool PhotometricMapping::checkDistance(const Transf & xi1,
                                         const Transf & xi2) const
{
    Transf delta = xi1.inverseCompose(xi2);
    double r = delta.rot().squaredNorm();
    double d = delta.trans().squaredNorm();
    return (r < _params.angularThreshSq and
            d < _params.distThreshSq);
}

bool PhotometricMapping::checkDistance(const Transf & xi) const
{
    double r = xi.rot().squaredNorm();
    double d = xi.trans().squaredNorm();
    return (r < _params.angularThreshSq and
            d < _params.distThreshSq);
}

Transf PhotometricMapping::localizeMI()
{
    ScalePhotometric localizer(5, _camera);
    localizer.setVerbosity(0);
    localizer.setXiBaseCam(_xiBaseCam);
    localizer.setTargetImage(_frameVec[_mapIdx].img);
    localizer.setBaseImage(_interFrame.img);
    localizer.setDepth(_depth);
    
    Transf & xiFr = _interFrame.xi;
    const Transf & xiMap = _frameVec[_mapIdx].xi;
    Transf xiFrMap = localizer.computePoseMI(xiFr.inverseCompose(xiMap));
    xiFr = xiMap.composeInverse(xiFrMap);
}

Transf PhotometricMapping::localizePhoto(const Mat8u & img)
{
    _localizer.setDepth(_depth);
    _localizer.setTargetImage(img);
    
    //estimated using only wheel odometry measurements
    Transf zetaPrior = _xiLocalOld.inverseCompose(_xiLocal);
    
    _xiLocal = _localizer.computePose( _xiLocal );
    Transf zeta = _xiLocalOld.inverseCompose(_xiLocal);
    double lengthPrior = zetaPrior.trans().norm();
    double length = zeta.trans().norm();
    double K = lengthPrior / length;
    if (length > 1e-2) //TODO put a meaningful threshold
    {
        
        if (K > 1.2) //too much divergence between WO and VO
        {
            cout << "WARNING : discard VO result, scale error" << endl;
            _xiLocal = _xiLocalOld.compose(zetaPrior);
        }
        else
        {
            zeta.trans() *= K;
            _xiLocal = _xiLocalOld.compose(zeta);
        }
    }
    _xiLocalOld = _xiLocal;
}








