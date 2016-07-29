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

#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_stereo.h"

using namespace std;
using namespace cv;

const int COLS = 1024;
const int ROWS = 768;
const double K = 0.3;

int main(int argc, char** argv)
{	
   
    array<double, 6> params = {0.5, 1, 250, 250, 512, 384};
    StereoParameters stereoParams;
    stereoParams.scale = 2;
    
    Transformation<double> T01(0.7, 0.1, 0.5, 0.1, -0.3, 0.5);
    Transformation<double> T0plane(0, 0, 1.5, 0, 0, 0);
    stereoParams.uMax = COLS;
    stereoParams.vMax = ROWS;
    
    EnhancedCamera camera(params.data());
    EnhancedStereo stereo(Transformation<double>(),
                &camera, &camera, stereoParams);
    
    
//     Init the localizer
    DepthMap depth0, depth1;
    stereo.generatePlane(T0plane, depth0,
         vector<Vector3d>{Vector3d(-0.5, -0.5, 0), Vector3d(0.5, -0.5, 0),
                          Vector3d(0.5, 0.5, 0), Vector3d(-0.5, 0.5, 0) } );
    
    stereo.generatePlane(T01.inverseCompose(T0plane), depth1,
         vector<Vector3d>{Vector3d(-0.5, -0.5, 0), Vector3d(0.5, -0.5, 0),
                          Vector3d(0.5, 0.5, 0), Vector3d(-0.5, 0.5, 0) } );
    
    DepthMap depth1wrap;
     
    DepthReprojector reprojector;
    
    Timer timer;
    reprojector.wrapDepth(depth0, depth1, T01, depth1wrap);
    cout << timer.elapsed() << endl;
    Mat32f img0(ROWS, COLS), img1(ROWS, COLS), img1wrap(ROWS, COLS);
    
    for (int y = 0; y < ROWS; y++)
    {
        for (int x = 0; x < COLS; x++)
        {
            img0(y, x) = depth0.nearest(x, y) * K;
            img1(y, x) = depth1.nearest(x, y) * K;
            img1wrap(y, x) = depth1wrap.nearest(x, y) * K;
        }
    }
    
    imshow("img0", img0);
    imshow("img1", img1);
    imshow("diff", img0 - img1wrap + 0.5);
    imshow("img1wrap", img1wrap);
    waitKey();
    
    return 0;
}



