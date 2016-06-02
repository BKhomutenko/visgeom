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

#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
#include <iomanip>

#include <boost/format.hpp>

#include "geometry.h"
#include "eucm.h"
#include "calibration/calibration.h"

using namespace cv;
using Eigen::JacobiSVD;
int main(int argc, char** argv) {

    vector<double> params1{0.5, 1, 500, 500, 500, 500};
    vector<double> params2{0.5, 1, 500, 500, 500, 500};

    Transformation<double> xi(0, 0, 0, 0, 0, 0);

    IntrinsicCameraCalibration<EnhancedProjector> calibLeft;
    IntrinsicCameraCalibration<EnhancedProjector> calibRight;


    ExtrinsicCameraCalibration<EnhancedProjector> calibStereo;

    if (calibStereo.initialize(argv[1], argv[2], argv[3]))
    {
        calibStereo.compute(params1, params2, xi);
        cout << "Extrinsic : " << endl;
        cout << xi << endl;
        cout << "Intrinsics left :" << endl;
        for (auto & p : params1)
        {
            cout << setw(12) << p;
        }
        cout << endl;
        cout << "Intrinsics right :" << endl;
        for (auto & p : params2)
        {
            cout << setw(12) << p;
        }
        cout << endl;
        calibStereo.residualAnalysis(params1, params2, xi);
        
    }

    return 0;

}
