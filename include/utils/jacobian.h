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
        Matrix23drm projJac;
        _camera->projectionJacobian(X2, projJac.data(), projJac.data() + 3);
        
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
        Matrix23drm projJac;
        _camera->projectionJacobian(X2, projJac.data(), projJac.data() + 3);
        
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

/*
computes the jacobian matrix dp / dxi23 if the transformation chain looks like
T12 T23
1 -- camera frame
2 -- intermediate frame
3 -- frame where the point is projected

inverted means that frame 2 follows frame 3, not precede

For example, if frames 1 and 3 coinside ,then xi13 will be Identity, and xi23 defines the camera pose
in a certain frame 2
*/
class InterJacobian
{
public:
    InterJacobian(const ICamera * camera,
            const Transformation<double> & xi13, const Transformation<double> & xi23, bool inverted) :
    _camera(camera->clone()),
    R12( xi13.rotMat() * xi23.rotMatInv() ),
    t13( xi13.trans() ),
    M12( R12 *interOmegaRot(xi23.rot()) ) 
    {
        if (inverted)
        {
            R12 *= -1; // to take into account the invertion of the kinematic screw
            M12 *= -1;
        }
        
    }
    
    // Point jacobian
    void dpdxi(const Vector3d & X1, double * dudxi, double * dvdxi)
    {
        Matrix23drm projJac;
        _camera->projectionJacobian(X1, projJac.data(), projJac.data() + 3);
        
        Vector3d t3X = X1 - t13; // projected into the first frame
        
        Map<Covector3d> dudtr(dudxi);
        Map<Covector3d> dudrot(dudxi + 3);
        dudtr = projJac.row(0) * R12;
        dudrot = -projJac.row(0) * hat(t3X) * M12;
        
        Map<Covector3d> dvdtr(dvdxi);
        Map<Covector3d> dvdrot(dvdxi + 3);
        dvdtr = projJac.row(1) * R12;
        dvdrot = -projJac.row(1) * hat(t3X) * M12;
    }
    
    //brightness jacobian
    void dfdxi(const Vector3d & X1, const Covector2d & grad, double * dfdxi)
    {
        Matrix23drm projJac;
        _camera->projectionJacobian(X1, projJac.data(), projJac.data() + 3);
        
        Covector3d dfdX = grad * projJac;
        
        Vector3d t3X = X1 - t13; // projected into the first frame
        
        Map<Covector3d> dfdtr(dfdxi);
        Map<Covector3d> dfdrot(dfdxi + 3);
        dfdtr = dfdX * R12;
        dfdrot = -dfdX * hat(t3X) * M12;
    }
    
private:
    ICamera * _camera;
    Matrix3d R12, M12;
    Vector3d t13;
    //R12 IS NOT NECESSARILY A ROTATION MATRIX
};



