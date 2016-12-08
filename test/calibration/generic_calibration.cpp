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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "geometry/geometry.h"
#include "projection/eucm.h"
#include "calibration/unified_calibration.h"

using boost::property_tree::ptree;
using boost::property_tree::read_json;

int main(int argc, char** argv) {
    
    GenericCameraCalibration calibration;
    
    for (int i = 1; i < argc; i++)
    {
        calibration.addResiduals(argv[i]);
    }
    
    calibration.compute();
    
    return 0;
}
