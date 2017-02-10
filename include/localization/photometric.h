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
Relative camera pose estimation based on photometric error and depth map
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"
#include "ceres.h"

#include "reconstruction/depth_map.h"
#include "localization/scale_space.h"
#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "localization/local_cost_functions.h"

//TODO add assertions ???
class ScalePhotometric
{
public:
    ScalePhotometric(int nScales, const ICamera * cam2) :
            scaleSpace1(nScales, true),
            scaleSpace2(nScales, false),
            camPtr2(cam2->clone()),
            _xiBaseCam(0, 0, 0, 0, 0, 0),
            verbosity(0) {}
            
           
    virtual ~ScalePhotometric()
    {
        delete camPtr2;
        camPtr2 = NULL;
    }
    
    void setXiBaseCam(const Transf & xiBaseCam) { _xiBaseCam = xiBaseCam; }
    void setNumberScales(int numScales)
    {
        scaleSpace1.setNumberScales(numScales);
        scaleSpace2.setNumberScales(numScales);
    }
    
    const DepthMap & depth() const { return depthMap; }
    DepthMap & depth() { return depthMap; }
    void setDepth(const DepthMap & newDepth) { depthMap = newDepth; }
    
    void computeBaseScaleSpace(const Mat32f & img1);
    
    void computePose(const Mat32f & img2, Transf & T12);
    
    void setVerbosity(int newVerbosity) { verbosity = newVerbosity; }
    
    void computePoseMI(const Mat32f & img2, Transf & T12);

    //TODO make enum for choosing the camera
    array<double, 6> covarianceEigenValues(const int scaleIdx, 
            const Transf T12, bool baseValues);

    //FIXME temporary function
    void wrapImage(const Mat32f & src, Mat32f & dst, const Transf T12) const;
private:
    // scaleSpace2 must be initialized
    void computePose(int scaleIdx, Transf & T12);
    void computePoseAuto(int scaleIdx, Transf & T12);
    void computePoseMI(int scaleIdx, Transf & T12);
    PhotometricPack initPhotometricData(int scaleIdx);

    Transf _xiBaseCam;

    BinaryScalSpace scaleSpace1;
    BinaryScalSpace scaleSpace2;
    ICamera * camPtr2;
    DepthMap depthMap;
    
    //TODO make a parameter structure
    // minimal squared norm of gradient for a pixel to be accepted
    const double GRAD_THRESH = 50;
    const double GRAD_MAX = 255;
    int verbosity;
};



