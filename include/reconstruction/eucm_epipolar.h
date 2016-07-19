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
Semi-global block matching algorithm for non-rectified images
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/curve_rasterizer.h"

class EnhancedEpipolar
{
public:
    EnhancedEpipolar(Transformation<double> T12,
        const double * params2, const int numberSteps, int verbosity = 0) :
        // initialize the members
        Transform12(T12),
        camera2(params2),
        step(4. / numberSteps),
        nSteps(numberSteps),
        verbosity(verbosity)
    {
        initialize();
    }
    
    const Polynomial2 & getCurve(Vector3d X) const { return epipolarVec[index(X)]; }
    
    void initialize();
    
private:
    
    void prepareCamera();
    void computePolynomial(Vector3d plane, Polynomial2 & surf) const;
    int index(Vector3d X) const;
    
    // variables for the epipolar computations
    // initialized with prepareCamera()
    double alpha;
    double beta;
    double fu;
    double fv;
    double u0;
    double v0;
    
    double gamma;
    double ag;
    double a2b;
    double fufv;
    double fufu;
    double fvfv;
    
    //must be computed before the epipolar computations
    Vector2d epipole;
    
    // pose of the first to the second camera
    Transformation<double> Transform12;  
    EnhancedCamera camera2;
   
    // total number of plains
    // must be even
    int nSteps;  
    
    // angular distance between two plains that represent the epipolar lines
    double step;
    
    // the basis in which the input vector is decomposed
    Vector3d xBase, yBase, zBase;
    
    // the epipolar curves represented by polynomial functions
    // epipolarVec[0] corresponds to base rotated about t by -pi/2
    std::vector<Polynomial2> epipolarVec;  
    
    int verbosity;
};

