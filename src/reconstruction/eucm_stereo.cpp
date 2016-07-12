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
#include "std.h"
#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"
#include "reconstruction/depth_map.h"

void EnhancedStereo::computeEpipole()
{
    Vector3d t21 = Transform12.transInv();
    cam2.projectPoint(t21, epipole);
    epipolePx[0] = round(epipole[0]);
    epipolePx[1] = round(epipole[1]);
}

CurveRasterizer<int, Polynomial2> EnhancedStereo::getCurveRasteriser(int idx)
{
    Vector2i pinfPx = pinfPxVec[idx];
    return CurveRasterizer<int, Polynomial2>(pinfPx, epipolePx, epipolarVec[idx]);
}

//TODO reconstruct the depth points, not everything
void EnhancedStereo::computeReconstructed()
{
    pointVec1.reserve(params.smallHeight*params.smallWidth);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            pointVec1.emplace_back(params.uBig(u), params.vBig(v));
        }
    }
    cam1.reconstructPointCloud(pointVec1, reconstVec);
}

// computes epipolarDirectionVec -- computed by shifting the reconstructed points in the direction 
// of motion infinitesimally and projecting them back
void EnhancedStereo::computeEpipolarDirections()
{
    epipolarDirectionVec.clear();
    Vector3d t = Transform12.trans();
    t = t.normalized() * 0.001;
    for (int i = 0; i < pointVec1.size(); i++)
    {
        Vector3d Xt = reconstVec[i] - t;
        Vector2d pt;
        cam1.projectPoint(Xt, pt);
        pt -= pointVec1[i];
        epipolarDirectionVec.push_back(pt.normalized());
    }
}

void EnhancedStereo::computeRotated()
{
    Transform12.inverseRotate(reconstVec, reconstRotVec);
}

void EnhancedStereo::computePinf()
{
    cam2.projectPointCloud(reconstRotVec, pinfVec);
    pinfPxVec.resize(pinfVec.size());
    for (int i = 0; i < pinfVec.size(); i++)
    {
        pinfPxVec[i][0] = round(pinfVec[i][0]);
        pinfPxVec[i][1] = round(pinfVec[i][1]);
    }
}

void EnhancedStereo::computeEpipolarCurves()
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeEpipolarCurves" << endl;
    Vector3d t21 = Transform12.transInv();
    epipolarVec.clear();
    const double & alpha = cam2.params[0];
    const double & beta = cam2.params[1];
    const double & fu = cam2.params[2];
    const double & fv = cam2.params[3];
    const double & u0 = cam2.params[4];
    const double & v0 = cam2.params[5];
    
    double gamma = 1 - alpha;
    double ag = (alpha - gamma);
    double a2b = alpha*alpha*beta;
    
    double fufv = fu * fv;
    double fufu = fu * fu;
    double fvfv = fv * fv;
    
    for (const auto & pt : reconstRotVec)
    {
        epipolarVec.emplace_back();
        Polynomial2 & surf = epipolarVec.back();
        Vector3d plane = pt.cross(t21);
        const double & A = plane[0];
        const double & B = plane[1];
        const double & C = plane[2];
        double AA = A * A;
        double BB = B * B;
        double CC = C * C;
        double CCfufv = CC * fufv;
        if (CCfufv/(AA + BB) < 0.5) // the curve passes through the projection center
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
}

void EnhancedStereo::traceEpipolarLine(int x, int y, Mat8u & out)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::traceEpipolarLine" << endl;
    int idx = getLinearIndex(x, y);
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(idx);
    
    int count = (pinfPxVec[idx] - epipolePx).norm();
    for (int i = 0; i < count; i++)
    {
        out(raster.y, raster.x) = 0;
        out(raster.y+1, raster.x) = 0;
        out(raster.y, raster.x+1) = 0;
        out(raster.y+1, raster.x+1) = 0;
        raster.step();
    }
}

void EnhancedStereo::createBuffer()
{
    if (params.verbosity > 1) cout << "EnhancedStereo::createBuffer" << endl;
    int bufferWidth = params.smallWidth*params.dispMax;
    if (errorBuffer.cols != bufferWidth or errorBuffer.rows != params.smallHeight)
    {
        errorBuffer = Mat8u(Size(bufferWidth, params.smallHeight));
    }
    if (tableauLeft.cols != bufferWidth or tableauLeft.rows != params.smallHeight)
    {
        tableauLeft = Mat32s(Size(bufferWidth, params.smallHeight));
    }
    if (tableauRight.cols != bufferWidth or tableauRight.rows != params.smallHeight)
    {
        tableauRight = Mat32s(Size(bufferWidth, params.smallHeight));
    }
    if (tableauTop.cols != bufferWidth or tableauTop.rows != params.smallHeight)
    {
        tableauTop = Mat32s(Size(bufferWidth, params.smallHeight));
    }
    if (tableauBottom.cols != bufferWidth or tableauBottom.rows != params.smallHeight)
    {
        tableauBottom = Mat32s(Size(bufferWidth, params.smallHeight));
    }
    if (smallDisparity.cols != params.smallWidth or smallDisparity.rows != params.smallHeight)
    {
        smallDisparity = Mat8u(Size(params.smallWidth, params.smallHeight));

    }
    if (params.verbosity > 2) 
    {
        cout << "    small disparity size: " << smallDisparity.size() << endl;
    }
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, 
        const Mat8u & img2,
        Mat8u & disparity)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    smallDisparity.copyTo(disparity);
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, 
        const Mat8u & img2,
        DepthMap & depth)
{
    computeCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    //TODO fuse the following code with computeDistance
    depth = DepthMap(&cam1, params.smallHeight, params.smallWidth,
            params.u0, params.v0, params.scale);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            depth.at(u, v) = computeDistance(u, v);
        }
    }
}

//TODO possibly change the step length along the curve to reduce the computation cost
void EnhancedStereo::computeCurveCost(const Mat8u & img1, const Mat8u & img2)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeCurveCost" << endl;
    
    const int HALF_LENGTH = max(params.scale - 1, 1);
    const int LENGTH = HALF_LENGTH * 2 + 1;
    
    // compute the weights for matching cost
    vector<int> weightVec(LENGTH);
    weightVec[HALF_LENGTH] = HALF_LENGTH + 1;
    for (int i = 0; i < HALF_LENGTH; i++)
    {
        weightVec[i] = i + 1;
        weightVec[LENGTH - i - 1] = i + 1;
    }
    int normalizer = accumulate(weightVec.begin(), weightVec.end(), 0);
    
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            int idx = getLinearIndex(u, v);
            
            // compute the local image descriptor -- a piece of the epipolar curve on the first image
            vector<uint8_t> descriptor;
            descriptor.reserve(LENGTH);
            Vector2d pt = pointVec1[idx];
            Vector2d dir = epipolarDirectionVec[idx];
            if (dir != dir)
            {
                uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;
                for (int d = 0; d < params.dispMax; d++, outPtr++)
                {
                    *outPtr = 0;
                }
                continue;
            }
            for (int i = -HALF_LENGTH; i <= HALF_LENGTH; i++)
            {
                Vector2d shifted = pt + i * dir;
                descriptor.push_back( bilinear(img1, shifted[0], shifted[1]) );
            }
            
            //sample img2 along the epipolar curve
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(idx);
            raster.steps(-HALF_LENGTH);
            vector<uint8_t> sampleVec(params.dispMax + LENGTH - 1, 0);
            for (int i = 0; i  < params.dispMax + LENGTH - 1; i++, raster.step())
            {
                if (raster.y < 0 or raster.y >= img2.rows 
                    or raster.x < 0 or raster.x >= img2.cols) sampleVec[i] = 0;
                else sampleVec[i] = img2(raster.y, raster.x);
            }
            
            //compute the bias;
            int sum1 = accumulate(descriptor.begin(), descriptor.end(), 0);
            
            // fill up the cost buffer
            uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;
            for (int d = 0; d < params.dispMax; d++, outPtr++)
            {
                int acc = 0;
                int sum2 = accumulate(sampleVec.begin() + d, sampleVec.begin() + d + LENGTH, 0);
                int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
                for (int i = 0; i < LENGTH; i++)
                {
                    acc += abs(descriptor[i] - sampleVec[d + i] + bias) * weightVec[i];
                }
                *outPtr = acc / normalizer;
            }
            
        }
    }
}
          
void EnhancedStereo::computeCost(const Mat8u & img1, const Mat8u & img2)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeCost" << endl;
    Mat8u img2remap(Size(params.scale + params.dispMax - 1, params.scale));
    Mat32s integral1, integral2;
    integral(img1, integral1);
    int scaleSquared = params.scale * params.scale;
    int hblock = int(params.scale) / 2;
    
    
    double radius = (params.scale - 1) / 2.;
    double centerShift;
    if (int(round(params.scale)) % 2) centerShift = 1;
    else centerShift = 0.5;
    
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            int idx = getLinearIndex(u, v);
            //FIXME this makes the algorithm dependent on the motion direction
            
            
//            Point pPxinf = pinfPxVec[idx];
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(idx);
            
            // Remap the band
            img2remap.setTo(0);
            
            // the right end
            int uBase = round(raster.x + centerShift);
            int vBase = round(raster.y - radius);
            int dstBase = params.dispMax + hblock;
            for (int i = 0; i < hblock; i++)
            {
                int u2 = uBase + i;
                if (u2 < 0 or u2 >= img2.cols) continue;
                for (int j = 0; j < params.scale; j++)
                {
                    int v2 = vBase + j;
                    if (v2 < 0 or v2 >= img2.rows) continue;
                    img2remap(j, dstBase + i) = img2(v2, u2);
                }
            }
            
            // the middle and the left
            for (int i = params.dispMax - 1 + hblock; i  >= 0; i--, raster.step())
            {
                int u2 = round(raster.x + centerShift - 1);
                int vBase = round(raster.y - radius);
                if (u2 < 0 or u2 >= img2.cols) continue;
                for (int j = 0; j < params.scale; j++)
                {
                    int v2 = vBase + j;
                    if (v2 < 0 or v2 >= img2.rows) continue;
                    img2remap(j, i) = img2(v2, u2);
                }
            }
                       
            
            //compute bias
            int u1 = params.uBig(u) - hblock;
            int v1 = params.vBig(v) - hblock;
            int bias1 = integral1(v1, u1) + integral1(v1 + params.scale, u1 + params.scale) -
                         integral1(v1 + params.scale, u1) - integral1(v1, u1 + params.scale);
//            
//            
            integral(img2remap, integral2);
            
            // compute the actual error
            uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;
            for (int i = params.dispMax - 1; i >= 0; i--, outPtr++)
            {
                int bias = integral2(params.scale, i + params.scale) - integral2(params.scale, i);
                bias = (bias - bias1) / scaleSquared;
                bias = min(10, max(-10, bias));
                int acc = 0;
                for (int x2 = 0; x2 < params.scale; x2++)
                {
                    for (int x1 = 0; x1 < params.scale; x1++)
                    {
                        acc += abs(img1(v1 + x2, u1 + x1) - 
                            img2remap(x2, i + x1) + bias);
                    }
                }
                
                *outPtr = acc / scaleSquared;
            }
        }
    }
} 

void EnhancedStereo::computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost)
{
    int bestCost = inCost[0];
    for (int i = 1; i < params.dispMax; i++)
    {
        bestCost = min(bestCost, inCost[i]);
    }
    int & val0 = outCost[0];
    val0 = inCost[0];
    val0 = min(val0, inCost[1] + params.lambdaStep);
    val0 = min(val0, bestCost + params.lambdaJump);
    val0 += error[0];
    for (int i = 1; i < params.dispMax-1; i++)
    {
        int & val = outCost[i];
        val = inCost[i];
        val = min(val, inCost[i + 1] + params.lambdaStep);
        val = min(val, inCost[i - 1] + params.lambdaStep);
        val = min(val, bestCost + params.lambdaJump);
        val += error[i];
    }
    int & vald = outCost[params.dispMax - 1];
    vald = inCost[params.dispMax - 1];
    vald = min(vald, inCost[params.dispMax - 2] + params.lambdaStep);
    vald = min(vald, bestCost + params.lambdaJump);
    vald += error[params.dispMax - 1];
}

void EnhancedStereo::computeDynamicProgramming()
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeDynamicProgramming" << endl;
    if (params.verbosity > 1) cout << "    left" << endl;
    // left tableau init
    for (int v = 0; v < params.smallHeight; v++)
    {
        int * tableauRow = (int *)(tableauLeft.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        // init the first row
        copy(errorRow, errorRow + params.dispMax, tableauRow);
        // fill up the tableau
        for (int u = 1; u < params.smallWidth; u++)
        {
            computeDynamicStep(tableauRow + (u - 1)*params.dispMax,
                    errorRow + u*params.dispMax, tableauRow + u*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    right" << endl;  
    // right tableau init
    for (int v = 0; v < params.smallHeight; v++)
    {
        int * tableauRow = (int *)(tableauRight.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        int base = (params.smallWidth - 1) * params.dispMax;
        copy(errorRow + base, errorRow + base + params.dispMax, tableauRow + base);
        
        for (int u = params.smallWidth - 2; u >= 0; u--)
        {
            computeDynamicStep(tableauRow + (u + 1)*params.dispMax, 
                    errorRow + u*params.dispMax, tableauRow + u*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    top" << endl;
    // top-down tableau init
    for (int u = 0; u < params.smallWidth; u++)
    {
        auto tableauCol = tableauTop(Rect(u*params.dispMax, 0, params.dispMax, params.smallHeight));
        auto errorCol = errorBuffer(Rect(u*params.dispMax, 0, params.dispMax, params.smallHeight));
        copy(errorCol.data, errorCol.data + params.dispMax, (int*)(tableauCol.data));
        for (int v = 1; v < params.smallHeight; v++)
        {
            computeDynamicStep((int*)(tableauCol.row(v-1).data), 
                    errorCol.row(v).data,
                    (int*)(tableauCol.row(v).data));
        }
    }
    if (params.verbosity > 1) cout << "    bottom" << endl;
    // bottom-up tableau init
    for (int u = 0; u < params.smallWidth; u++)
    {
        auto tableauCol = tableauBottom(Rect(u*params.dispMax, 0, params.dispMax, params.smallHeight));
        auto errorCol = errorBuffer(Rect(u*params.dispMax, 0, params.dispMax, params.smallHeight));
        int vLast = params.smallHeight - 1;
        copy(errorCol.row(vLast).data, 
                errorCol.row(vLast).data + params.dispMax, 
                (int*)(tableauCol.row(vLast).data));
        for (int v = params.smallHeight - 2; v >= 0; v--)
        {
            computeDynamicStep((int*)(tableauCol.row(v+1).data), 
                    errorCol.row(v).data,
                    (int*)(tableauCol.row(v).data));
        }
    }
    
}

void EnhancedStereo::reconstructDisparity()
{
    if (params.verbosity > 0) cout << "EnhancedStereo::reconstructDisparity" << endl;
    Mat16s errFinalMat(smallDisparity.size());
    for (int v = 0; v < params.smallHeight; v++)
    {
        int* dynRow1 = (int*)(tableauLeft.row(v).data);
        int* dynRow2 = (int*)(tableauRight.row(v).data);
        int* dynRow3 = (int*)(tableauTop.row(v).data);
        int* dynRow4 = (int*)(tableauBottom.row(v).data);
        uint8_t* errRow = (errorBuffer.row(v).data);
        for (int u = 0; u < params.smallWidth; u++)
        {
            int bestCost = 100000;
            uint8_t & bestDisp = smallDisparity(v, u);
            int16_t & bestErr = errFinalMat(v, u);
            bestDisp = 0;
            for (int d = 0; d < params.dispMax; d++)
            {
                int base = u * params.dispMax;
//                if (errRow[base + d] > 5) continue;
//                int acc = errRow[base + d];
                int acc =dynRow1[base + d] + dynRow2[base + d] 
                        + dynRow3[base + d] + dynRow4[base + d] - 2*errRow[base + d];
                if (bestCost > acc)
                {
                    bestCost = acc;
                    bestDisp = d;
                    bestErr = acc;
                }
            }
            if (params.verbosity > 3) cout << "    best error: " << int(bestErr) << endl;
        }
        
    }
}

//TODO make error codes
bool EnhancedStereo::triangulate(double x1, double y1, double x2, double y2, Vector3d & X)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    //Vector3d v1n = v1 / v1.norm(), v2n = v2 / v2.norm();
    Vector3d v1, v2;
    if (not cam1.reconstructPoint(Vector2d(x1, y1), v1) or 
        not cam2.reconstructPoint(Vector2d(x2, y2), v2) )
    {
        if (params.verbosity > 2) 
        {
            cout << "    not reconstructed " << Vector2d(x1, y1).transpose(); 
            cout << " # " << Vector2d(x2, y2).transpose() << endl;
        }
        X = Vector3d(0, 0, 0);
        return false;
    }
    Vector3d t = Transform12.trans();
    v2 = Transform12.rotMat() * v2;
    if (params.verbosity > 3) 
    {
        cout << "    pt1: " << x1 << " " << y1 << endl;
        cout << "    x1: " << v1.transpose() << endl;
        cout << "    pt2: " << x2 << " " << y2 << endl;
        cout << "    x2: " << v2.transpose() << endl;
    }
    double v1v2 = v1.dot(v2);
    double v1v1 = v1.dot(v1);
    double v2v2 = v2.dot(v2);
    double tv1 = t.dot(v1);
    double tv2 = t.dot(v2);
    double delta = -v1v1 * v2v2 + v1v2 * v1v2;
    if (abs(delta) < 1e-10) // TODO the constant to be revised
    {
        if (params.verbosity > 2) 
        {
            cout << "    not triangulated " << abs(delta) << " " << (abs(delta) < 1e-10) << endl;
        }
        X = Vector3d(0, 0, 0);
        return false;
    }
    double l1 = (-tv1 * v2v2 + tv2 * v1v2)/delta;
    double l2 = (tv2 * v1v1 - tv1 * v1v2)/delta;
    X = (v1*l1 + t + v2*l2)*0.5;
    return true;
}

void EnhancedStereo::computeDistance(Mat32f & distance)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeDistance(Mat32f &)" << endl;
    distance.create(params.smallHeight, params.smallWidth);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            distance(v, u) = computeDistance(u, v);
        }
    }
}

double EnhancedStereo::computeDistance(int u, int v)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::computeDistance(int *)" << endl;
    
    int idx = getLinearIndex(u, v);
    int disparity = smallDisparity(v, u);
    if (disparity <= 0) 
    {
        return params.maxDistance;
    }
    
    // to compute point on the second image
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(idx);
    raster.steps(disparity);
    
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    Vector3d X;
    if (triangulate(pt1[0], pt1[1], raster.x, raster.y, X))
    {
        return X.norm();
    }
    else 
    {
        return 0;
    }
}

void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        Mat32f & distance, const Vector3dVec & polygonVec)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::generatePlane" << endl;
    distance.create(params.smallHeight, params.smallWidth);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            distance(v, u) = 0;
            Vector3d vec; // the direction vector
            if (not cam1.reconstructPoint(Vector2d(params.uBig(u), params.vBig(v)), vec)) continue;
            double zvec = z.dot(vec);
            if (zvec < 1e-3) 
            {
                continue;
            }
            bool inside = true;
            for (int i = 0; i < polygonCamVec.size(); i++)
            {
                int j = (i + 1) % polygonCamVec.size();
                Vector3d normal = polygonCamVec[i].cross(polygonCamVec[j]);
                if (vec.dot(normal) < 0)
                {
                    inside = false;
                    break;
                }
            }
            if (not inside) continue;
            double tz = t.dot(z);
            double alpha = tz / zvec;
            vec *= alpha;
            distance(v, u) = vec.norm();
        }
    }
}

void EnhancedStereo::upsampleDisparity(const Mat8u & img1, Mat8u & disparity)
{
    cout << smallDisparity.size() << endl;
    smallDisparity.copyTo(disparity);
//    resize(smallDisparity, disparity, Size(0, 0), params.scale, params.scale, 0);
}

