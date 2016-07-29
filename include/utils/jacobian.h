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

#pragma once

#include "eigen.h"

#include "camera/generic_camera.h"
#include "geometry/geometry.h"

// A class for projection jacobian computation
// dp/dt = J dxi/dt

//TODO rewrite the explanation 
/*
Jacobian matrix of wrapped image brightness with respect
to xi -- the transformation from camera frame 2 to camera 
frame 1 transformation.
Transformation is parametrized using uTheta angle-axis

J = grad(img) * dp/dX * [ -I  hat(X) ] * L_uTheta 

- grad(img) is computed with the interpolation

- dp/dX is provided by the camera

- [ -I  hat(X) ] represents the relation between motion of a spatial point 
and camera's kinematic screw. X and [  v  ] must be expressed in frame 2
                                    [omega]
dX/dt = [ -I  hat(X) ] * [  v  ]
                         [omega]

- L_uTheta is a mapping between dxi/dt and [  v  ] = V in frame 2
                                           [omega]     
L_uTheta =  [ R21      0        ]
            [  0   R21*M_uTheta ]
            
V = L_uTheta * dxi/dt
*/
class CameraJacobian
{
public:
    CameraJacobian(const ICamera * camera, const Transformation<double> & T12) :
        _camera(camera->clone()),
        R21( T12.rotMatInv() ),
        M21( R21 * interOmegaRot(T12.rot()) ) {}
    
    // Point jacobian
    void dpdxi(const Vector3d & X2, double * dudxi, double * dvdxi)
    {
        Matrix23d projJac;
        _camera->projectionJacobian(X2, projJac);
        
        Map<Covector3d> dudtr(dudxi);
        Map<Covector3d> dudrot(dudxi + 3);
        dudtr = -projJac.row(0) * R21;
        dudrot = projJac.row(0) * hat(X2) * M21;
        
        Map<Covector3d> dvdtr(dvdxi);
        Map<Covector3d> dvdrot(dvdxi + 3);
        dvdtr = -projJac.row(1) * R21;
        dvdrot = projJac.row(1) * hat(X2) * M21;
    }
    
    //brightness jacobian
    void dfdxi(const Vector3d & X2, const Covector2d & grad, double * dfdxi)
    {
        Matrix23d projJac;
        _camera->projectionJacobian(X2, projJac);
        
        Covector3d dfdX = grad * projJac;
        
        Map<Covector3d> dfdtr(dfdxi);
        Map<Covector3d> dfdrot(dfdxi + 3);
        dfdtr = -dfdX * R21;
        dfdrot = dfdX * hat(X2) * M21;
    }
    
private:
    ICamera * _camera;
    Matrix3d R21, M21;
};

