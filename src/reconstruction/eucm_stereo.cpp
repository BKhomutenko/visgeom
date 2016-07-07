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
    epipolePx.x = round(epipole[0]);
    epipolePx.y = round(epipole[1]);
}

//TODO reconstruct the depth points, not everything
void EnhancedStereo::computeReconstructed()
{
    Vector2dVec imagePointVec;
    imagePointVec.reserve(params.smallHeight*params.smallWidth);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            imagePointVec.emplace_back(params.uBig(u), params.vBig(v));
        }
    }
    cam1.reconstructPointCloud(imagePointVec, reconstVec);
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
        pinfPxVec[i].x = round(pinfVec[i][0]);
        pinfPxVec[i].y = round(pinfVec[i][1]);
    }
        
}

void EnhancedStereo::computeEpipolarCurves()
{
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

void EnhancedStereo::traceEpipolarLine(Point pt, Mat8u & out)
{
    int idx = pt.y * cam1.width + pt.x;
    cout << idx << endl;
    cout << pinfPxVec.size() << endl;
    Point pinf = pinfPxVec[idx];
    cout << pinf << endl;
    CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
    
    int count = abs(pinf.x - epipolePx.x);
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
//    cout << "create" << endl;
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
        cout << smallDisparity.size() << endl;
    }
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, 
        const Mat8u & img2,
        Mat8u & disparity)
{
    computeCost(img1, img2);
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
            if (smallDisparity(v, u) == 0) 
            {
                depth.at(u, v) = 100;
                continue;
            }
            int idx = getLinearIdx(params.vBig(v), params.uBig(u));
            Point pinf = pinfPxVec[idx];
            CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
            raster.step(smallDisparity(v, u));
            Vector3d X = triangulate(params.uBig(u), params.vBig(v), raster.x, raster.y);
            depth.at(u, v) = X.norm();
        }
    }
}
          
void EnhancedStereo::computeCost(const Mat8u & img1, const Mat8u & img2)
{
    double T1 = 0, T2 = 0; // time profiling
//    cout << "cost" << endl;
    Mat8u img2remap(Size(params.scale + params.dispMax - 1, params.scale));
    Mat32s integral1, integral2;
    integral(img1, integral1);
    int scaleSquared = params.scale * params.scale;
    int hblock = int((params.scale - 1.) / 2.);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            int idx = getLinearIdx(u, v);
            //FIXME this makes the algorithm dependent on the motion direction
            uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;
            
//            Point pPxinf = pinfPxVec[idx];
            Vector2d pinf = pinfVec[idx];
            CurveRasterizer<Polynomial2> raster(pinf[0], pinf[1], epipole[0], epipole[1], epipolarVec[idx]);
            
            // Remap the band
            img2remap.setTo(0);
            
            // the right end
            int uBase = ceil(pinf[0]);
            int vBase = int(floor(pinf[1])) - hblock;
            int dstBase = params.dispMax - 1 + hblock;
            for (int i = 0; i <= hblock; i++)
            {
                int u2 = uBase + i;
                if (u2 < 0 or u2 >= img2.cols) continue;
                for (int j = 0; j <= params.scale; j++)
                {
                    int v2 = vBase + j;
                    if (v2 < 0 or v2 >= img2.rows) continue;
                    img2remap(j, dstBase + i) = img2(v2, u2);
                }
            }
            
            // the middle and the left
            for (int i = params.dispMax - 2 + hblock; i  >= 0; i--, raster.step())
            {
                int u2 = ceil(raster.x);
                int vBase = int(floor(raster.y)) - hblock;
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
            
            for (int i = 0; i < params.dispMax; i++, outPtr++)
            {
                int bias = integral2(params.scale, i + params.scale) - integral2(params.scale, i);
                bias = (bias - bias1) / scaleSquared;
                bias = min(10, max(-10, bias));
//                cout << bias << " ";
                int acc = 0;
                for (int x2 = -hblock; x2 <= hblock; x2++)
                {
                    for (int x1 = -hblock; x1 <= hblock; x1++)
                    {
                        acc += abs(img1(v1 + x2, u1- x1) - 
                            img2remap(hblock + x2, i + hblock + x1) - bias);
                    }
                }
                
                *outPtr = acc / scaleSquared;
            }
//            cout << endl;
        }
    }
//    cout << "read " << T1 / CLOCKS_PER_SEC << endl;
//    cout << "write " << T2 / CLOCKS_PER_SEC << endl;
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
//    cout << "left" << endl;
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
//    cout << "right" << endl;    
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
//    cout << "top" << endl;
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
//    cout << "bottom" << endl;
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
//            cout << int(bestErr) << endl;
        }
        
    }
    imshow("errFinal", errFinalMat);
//    medianBlur(smallDisparity, smallDisparity, 3);
}

Vector3d EnhancedStereo::triangulate(int x1, int y1, int x2, int y2)
{
    //Vector3d v1n = v1 / v1.norm(), v2n = v2 / v2.norm();
    Vector3d v1, v2;
    cam1.reconstructPoint(Vector2d(x1, y1), v1);
    cam2.reconstructPoint(Vector2d(x2, y2), v2);
    Vector3d t = Transform12.trans();
    v2 = Transform12.rotMat() * v2;
    double v1v2 = v1.dot(v2);
    double v1v1 = v1.dot(v1);
    double v2v2 = v2.dot(v2);
    double tv1 = t.dot(v1);
    double tv2 = t.dot(v2);
    double delta = -v1v1 * v2v2 + v1v2 * v1v2;
    if (abs(delta) < 1e-4) // TODO the constant to be revised
    {
        return Vector3d(0, 0, 0);
    }
    double l1 = (-tv1 * v2v2 + tv2 * v1v2)/delta;
    double l2 = (tv2 * v1v1 - tv1 * v1v2)/delta;
    return (v1*l1 + t + v2*l2)*0.5;
}

void EnhancedStereo::computeDistance(Mat32f & distance)
{
    distance.create(params.smallHeight, params.smallWidth);
    for (int v = 0; v < params.smallHeight; v++)
    {
        for (int u = 0; u < params.smallWidth; u++)
        {
            if (smallDisparity(v, u) == 0) 
            {
                distance(v, u) = 100;
                continue;
            }
            int idx = getLinearIdx(params.vBig(v), params.uBig(u));
            Point pinf = pinfPxVec[idx];
            CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
            raster.step(smallDisparity(v, u));
            Vector3d X = triangulate(params.uBig(u), params.vBig(v), raster.x, raster.y);
            distance(v, u) = X.norm();
        }
    }
}

void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        Mat32f & distance, const Vector3dVec & polygonVec)
{
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

