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
#include "projection/eucm.h"
#include "utils/curve_rasterizer.h"


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
    
    // prepare different directions
    Vector3dVec directionVec;
    for (int idx = 0; idx < nSteps; idx++)
    {
        if (idx < nSteps/2)
        {
            //tangent part
            double s = step * idx - 1;
            directionVec.emplace_back(xBase + s*yBase);
        }
        else
        {
            //cotangent part
            double c = step * (-idx + nSteps/2) + 1;
            directionVec.emplace_back(c * xBase + yBase);
        }
    }
    
    // ## camera 1 ##
    // prepare reused variables
    prepareCamera(CAMERA_1);
    epipolar1Vec.clear();
    epipolar1Vec.reserve(directionVec.size());
    for (const auto & dir : directionVec)
    {
        Vector3d plane = dir.cross(zBase);
        epipolar1Vec.emplace_back(computePolynomial(plane));
    }
    epipolar1Vec.emplace_back(epipolar1Vec.front());
    
    // ## camera 2 ##
    
    // the transformation to frame 2
    Matrix3d R21 = Transform12.rotMatInv();
    Vector3d t21n = R21 * zBase;
    
    // prepare reused variables
    prepareCamera(CAMERA_2);
    epipolar2Vec.clear();
    epipolar2Vec.reserve(directionVec.size());
    for (const auto & dir : directionVec)
    {
        Vector3d dir2 = R21 * dir;
        Vector3d plane = dir2.cross(t21n);
        epipolar2Vec.emplace_back(computePolynomial(plane));
    }
    epipolar2Vec.emplace_back(epipolar2Vec.front());
    
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

Polynomial2 EnhancedEpipolar::computePolynomial(Vector3d plane) const
{
    Polynomial2 surf;
    const double & A = plane[0];
    const double & B = plane[1];
    const double & C = plane[2];
    const double AA = A * A;
    const double BB = B * B;
    const double CC = C * C;
    const double CCfufv = CC * fufv;
    const double dd = CCfufv / (AA + BB); //squared distance to the projection center in pixels
    if ((AA + BB) > 0 and dd < 1.) // the curve passes through the projection center
    {
        
        surf.kuu = surf.kuv = surf.kvv = 0;
        surf.ku = A / fu;
        surf.kv = B / fv;
        
        const double normABinv = 1./sqrt(AA + BB);
        const double Cnorm = C / sqrt(AA + BB + CC);
        const double du = -A * Cnorm * normABinv * fu;
        const double dv = -B * Cnorm * normABinv * fv;
        surf.k1 = -(u0  + du)* A / fu - (v0 + dv) * B / fv;
        
//        surf.k1 = -u0 * A/ fu - v0 * B / fv;
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
    return surf;
}

void EnhancedEpipolar::prepareCamera(CameraIdx camIdx)
{
    const double * params = NULL;
    if (camIdx == CAMERA_1)
    {
        params = camera1->getParams();
    } 
    else if (camIdx == CAMERA_2) 
    {
        params = camera2->getParams();
    }
    epipole = epipoles->get(camIdx);
    activeCamera = camIdx;
    alpha = params[0];
    beta = params[1];
    fu = params[2];
    fv = params[3];
    u0 = params[4];
    v0 = params[5];
    
    // intermediate variables
    gamma = 1 - alpha;
    ag = (alpha - gamma);
    a2b = alpha*alpha*beta;
    fufv = fu * fv;
    fufu = fu * fu;
    fvfv = fv * fv;
}

void EnhancedEpipolar::traceEpipolarLine(int u, int v, Mat & out, CameraIdx camIdx, int count) const
{
    if (verbosity > 0) cout << "EnhancedStereo::traceEpipolarLine" << endl;
    
    CurveRasterizer<int, Polynomial2> * raster1 = NULL;
    Vector3d X1, X2;
    Vector2d pt;
    CameraIdx targetIdx = camIdx == CAMERA_1 ? CAMERA_2 : CAMERA_1;
    if (camIdx == CAMERA_1)
    { 
        if (not camera1->reconstructPoint(Vector2d(u, v), X1)) return;
        X2 = Transform12.rotMatInv() * X1;
        if (not camera2->projectPoint(X2, pt)) return;
    }
    else if (camIdx == CAMERA_2)
    {
        if (not camera2->reconstructPoint(Vector2d(u, v), X2)) return;
        X1 = Transform12.rotMat() * X2;
        if (not camera1->projectPoint(X1, pt)) return;
    }
    
    
    Vector2i pti = round(pt); //FIXME make a function. where to put?
    auto useInverted = epipoles->chooseEpipole(targetIdx, pti);
    Vector2i goal = epipoles->getPx(targetIdx, useInverted);
    raster1 = new CurveRasterizer<int, Polynomial2>(pti, goal, get(targetIdx, X1));
    if (useInverted & EPIPOLE_INVERTED) raster1->setStep(-1);
        
    CurveRasterizer<int, Polynomial2> * raster2 = new CurveRasterizer<int, Polynomial2>(*raster1);
    cv::circle(out, Point(raster1->u, raster1->v), 1, Scalar(128), -1);
    for (int i = 0; i < count; i++)
    {
        cv::circle(out, Point(raster1->u, raster1->v), 0, Scalar(128), -1);
//       cv::circle(out, Point(raster2->u, raster2->v), 0, Scalar(128), -1);
        raster1->step();
//        raster2->unstep();
    }
    delete raster1;
    delete raster2;
}
   
