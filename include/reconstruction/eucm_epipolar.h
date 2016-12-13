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
A class that computes the epipolar curve equations for a calibrated stereo system
*/

#pragma once

#include "ocv.h"
#include "std.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "projection/eucm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/epipoles.h"
#include "reconstruction/stereo_misc.h"

class EnhancedEpipolar
{
public:
   
    EnhancedEpipolar(const EnhancedCamera * cam1, const EnhancedCamera * cam2,
            const int numberSteps, int verbosity = 0) :
        // initialize the members
        camera1(cam1->clone()),
        camera2(cam2->clone()),
        step(4. / numberSteps),
        nSteps(numberSteps),
        verbosity(verbosity)
    {
    
            cout << "Desc constructor" << endl;
    }
    
    EnhancedEpipolar(const EnhancedCamera * cam1, const EnhancedCamera * cam2,
            const Transf & T12, const int numberSteps, int verbosity = 0) :
        // initialize the members
        Transform12(T12),
        camera1(cam1->clone()),
        camera2(cam2->clone()),
        epipoles(cam1, cam2, T12),
        step(4. / numberSteps),
        nSteps(numberSteps),
        verbosity(verbosity)
    {
        assert(T12.trans().squaredNorm() > 1e-10);
        initialize();
    }
    
    virtual ~EnhancedEpipolar()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
    }
    
    void setTransformation(const Transf & transf)
    {
        assert(transf.trans().squaredNorm() > 1e-10);
        Transform12 = transf;
        epipoles = StereoEpipoles(camera1, camera2, transf);
        initialize();
    }
    
    const Polynomial2 & getFirst(Vector3d X) const { return epipolar1Vec[index(X)]; }
    
    const Polynomial2 & getSecond(Vector3d X) const { return epipolar2Vec[index(X)]; }
    
    // draws an epipolar line  on the right image that corresponds to (x, y) on the left image
    void traceEpipolarLine(int u, int v, Mat & out, CameraIdx camIdx, int count = 150) const;
    
    void initialize();
    
    const StereoEpipoles & getEpipoles() const { return epipoles; }
    
private:
    
    void prepareCamera(CameraIdx camIdx);
    Polynomial2 computePolynomial(Vector3d plane) const;
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
    
    //stores the epipoles
    StereoEpipoles epipoles;
    
    //defines which camera is active right now
    CameraIdx activeCamera;
    
    Vector2d epipole;
    
    // pose of the first to the second camera
    Transf Transform12;  
    EnhancedCamera *camera1, *camera2;
   
    // total number of plains
    // must be even
    int nSteps;  
    
    // angular distance between two plains that represent the epipolar lines
    double step;
    
    // the basis in which the input vector is decomposed
    Vector3d xBase, yBase, zBase;
    
    // the epipolar curves represented by polynomial functions
    // for both cameras
    // epipolar1Vec[0] corresponds to base rotated about t by -pi/2
    std::vector<Polynomial2> epipolar2Vec;
    std::vector<Polynomial2> epipolar1Vec;
    
    int verbosity;
};

