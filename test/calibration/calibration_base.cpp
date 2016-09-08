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

#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "calibration/calibration.h"

int main(int argc, char** argv) {

//    vector<double> params{0.5, 1, 500, 500, 500, 500};
    vector<double> params{0.5,    1.,    500,     500,    500,    500};
    Transformation<double> TbaseCam, TorigGrid;

    BaseTransformationCalibration<EnhancedProjector> calibRobot;

    if (true or calibRobot.initialize(argv[1]))
    {
        calibRobot.compute(TbaseCam, TorigGrid);
        cout << "TbaseCam : " << endl;
        cout << TbaseCam << endl;
        cout << "TorigGrid :" << endl;
        cout << TorigGrid << endl;
        calibRobot.residualAnalysis(TbaseCam, TorigGrid);
       
        
//        TbaseCam = Transformation<double>(0, 0, 0, -1.5927, -0.0164123, 0.00799701);
//        auto R = TbaseCam.rotMat();
//        double thx = atan2(R(2, 1), R(2, 2));
//        double sthx = sin(thx);
//        double thy = atan2(-R(2, 0), R(2, 1) / sthx);
//       
//        cout << thx << " " << thy << endl;
//        Transformation<double> RotX(0, 0, 0, thx, 0, 0);
//        Transformation<double> RotY(0, 0, 0, 0, thy, 0);
//        cout << RotY.compose(RotX) << endl;
//        cout << TbaseCam.compose(RotY.compose(RotX).inverse()) << endl;
    }

    return 0;

}
