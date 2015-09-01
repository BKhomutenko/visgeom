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
Stereo vision definition.
*/

#ifndef _SPCMAP_VISION_H_
#define _SPCMAP_VISION_H_

//STL
#include <vector>
#include <algorithm>

//Eigen
#include <Eigen/Eigen>

#include "geometry.h"
#include "camera.h"

using namespace std;

enum CameraID {LEFT, RIGHT};

class StereoSystem
{
public:
    void projectPointCloud(const vector<Eigen::Vector3d> & src,
            vector<Eigen::Vector2d> & dst1, vector<Eigen::Vector2d> & dst2) const;

    void reconstructPointCloud(const vector<Eigen::Vector2d> & src1, const vector<Eigen::Vector2d> & src2,
            vector<Eigen::Vector3d> & dst) const;

    //TODO make smart constructor with calibration data passed
    StereoSystem(Transformation<double> & p1, Transformation<double> & p2,
            ICamera & c1, ICamera & c2)
            : TbaseCam1(p1), TbaseCam2(p2), cam1(c1.clone()), cam2(c2.clone()) {}

    ~StereoSystem();

    static bool triangulate(const Eigen::Vector3d & v1, const Eigen::Vector3d & v2,
            const Eigen::Vector3d & t,  Eigen::Vector3d & X);
    void reconstruct2(const Eigen::Vector2d & p1,
            const Eigen::Vector2d & p2,
            Eigen::Vector3d & X) const;
    Transformation<double> TbaseCam1;  // pose of the left camera in the base frame
    Transformation<double> TbaseCam2;  // pose of the right camera in the base frame
    ICamera * cam1, * cam2;
};

void computeEssentialMatrix(const vector<Eigen::Vector3d> & xVec1,
        const vector<Eigen::Vector3d> & xVec2,
        Eigen::Matrix3d & E);

#endif
