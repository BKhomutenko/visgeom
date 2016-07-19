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
#include "reconstruction/eucm_epipolar.h"

#include "timer.h"
#include "io.h"
#include "std.h"
#include "eigen.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/curve_rasterizer.h"


void EnhancedEpipolar::initialize()        
{
    assert(nSteps % 2 == 0);
    if (verbosity > 0) cout << "EnhancedEpipolar::initialize" << endl;
    Timer timer;
    // compute the epipolar basis
    // it is used to quickly access epipolar lines by a direction vector
    
    // translation vector defines z
    zBase = -Transform12.trans().normalized();
    
    // find a vector perpendicular to z
    if (zBase[2]*zBase[2] > zBase[0]*zBase[0] + zBase[1]*zBase[1])
    {
        xBase << 1, 0, 0;
    }
    else
    {
        xBase << 0, 0, 1;
    }
    Matrix3d orthProjector = Matrix3d::Identity() - zBase*zBase.transpose();
    xBase = orthProjector * xBase;
    xBase.normalize();
    
    // the last component of the basis
    yBase = zBase.cross(xBase);
    
    // the 
    Matrix3d R21 = Transform12.rotMatInv();
    Vector3d t21n = R21 * zBase;
    
    // prepare reused variables
    camera2.projectPoint(t21n, epipole);
    prepareCamera();
    
    for (int idx = 0; idx < nSteps; idx++)
    {
        Vector3d X;
        if (idx < nSteps/2)
        {
            //tangent part
            double s = step * idx - 1;
            X = xBase + s*yBase;
        }
        else
        {
            //cotangent part
            double c = step * (-idx + nSteps/2) + 1;
            X = c * xBase + yBase;
        }
        X = R21 * X;
        epipolarVec.emplace_back();
        Vector3d plane = X.cross(t21n);
        computePolynomial(plane, epipolarVec.back());
    }
    epipolarVec.emplace_back(epipolarVec.front());
    if (verbosity > 1) cout << "    epipolar init time : " << timer.elapsed() << endl;
}


int EnhancedEpipolar::index(Vector3d X) const
{
    double c = X.dot(xBase);
    double ac = abs(c);
    double s = X.dot(yBase);
    double as = abs(s);
    if (ac + as < 1e-4) //TODO check the constant 
    {
        return 0;
    }
    else if (ac > as)
    {
        return round((s/c + 1) / step);
    }
    else
    {
        return round((1 - c/s) / step) + nSteps/2;
    }
}

void EnhancedEpipolar::computePolynomial(Vector3d plane, Polynomial2 & surf) const
{
    const double & A = plane[0];
    const double & B = plane[1];
    const double & C = plane[2];
    double AA = A * A;
    double BB = B * B;
    double CC = C * C;
    double CCfufv = CC * fufv;
    if (CCfufv/(AA + BB) < 2e-1) // the curve passes through the projection center
    {
        surf.kuu = surf.kuv = surf.kvv = 0;
        surf.ku = A/fu;
        surf.kv = B/fv;
        surf.k1 = -u0*A/fu - v0*B/fv;
    }
    else
    {
        // compute first 4 coefficients directly
        surf.kuu = (AA*ag + CC*a2b)/(CC*fufu);  // kuu
        surf.kuv = 2*A*B*ag/(CCfufv);  // kuv
        surf.kvv = (BB*ag + CC*a2b)/(CC*fvfv);  // kvv
        surf.ku = 2*(-(AA*fv*u0 + A*B*fu*v0)*ag - 
                        A*C*fufv*gamma - CC*a2b*fv*u0)/(CCfufv*fu);  // kv
        surf.kv = 2*(-(BB*fu*v0 + A*B*fv*u0)*ag - 
                        B*C*fufv*gamma - CC*a2b*fu*v0)/(CCfufv*fv);  // kv
                        
        // the last one is computed using the fact that
        // the epipolar curves pass through the epipole
        surf.k1 = -(surf.kuu*epipole[0]*epipole[0] 
                        + surf.kuv*epipole[0]*epipole[1] 
                        + surf.kvv*epipole[1]*epipole[1] 
                        + surf.ku*epipole[0] + surf.kv*epipole[1]);
    }
}

void EnhancedEpipolar::prepareCamera()
{
    alpha = camera2.params[0];
    beta = camera2.params[1];
    fu = camera2.params[2];
    fv = camera2.params[3];
    u0 = camera2.params[4];
    v0 = camera2.params[5];
    
    // intermediate variables
    gamma = 1 - alpha;
    ag = (alpha - gamma);
    a2b = alpha*alpha*beta;
    fufv = fu * fv;
    fufu = fu * fu;
    fvfv = fv * fv;
}
   
