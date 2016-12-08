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

#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"

#include "projection/eucm.h"
#include "reconstruction/depth_map.h"

using namespace std;
using namespace cv;

const int COLS = 1024;
const int ROWS = 768;
const double K = 0.3;

int main(int argc, char** argv)
{	
   
    array<double, 6> params = {0.5, 1, 250, 250, 512, 384};
    ScaleParameters scaleParams;
    scaleParams.scale = 2;
    
    Transformation<double> T01(0.7, 0.1, 0.5, 0.1, -0.3, 0.5);
    Transformation<double> T0plane(0, 0, 0.5, 0, 0, 0);
    scaleParams.uMax = COLS;
    scaleParams.vMax = ROWS;
    scaleParams.setEqualMargin();
    EnhancedCamera camera(params.data());
    
//     Init the localizer
    DepthMap depth = DepthMap::generatePlane(&camera, scaleParams, T0plane,
         vector<Vector3d>{Vector3d(-0.5, -0.5, 0), Vector3d(0.5, -0.5, 0),
                          Vector3d(0.5, 0.5, 0), Vector3d(-0.5, 0.5, 0) } );
    
    Timer timer;
    MHPack pack;
    depth.reconstruct(pack, SIGMA_VALUE);
    DepthMap depthProjected(&camera, scaleParams);
    depthProjected.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
    for (int i = 0; i < pack.cloud.size(); i++)
    {
        depthProjected.pushHypothesis(pack.cloud[i], pack.sigmaVec[i]);
    }
    cout << timer.elapsed() << endl;

    Mat32f img0, img1;
    
    depth.toMat(img0);
    depthProjected.toMat(img1);
    
    imshow("img0", img0);
    imshow("img1", img1);
    imshow("diff", img0 - img1 + 0.5);
    waitKey();
    
    return 0;
}



