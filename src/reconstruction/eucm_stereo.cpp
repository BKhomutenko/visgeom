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
#include "filter.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"
#include "reconstruction/depth_map.h"

void EnhancedStereo::computeEpipole()
{
    if (not camera1->projectPoint(Transform12.trans(), epipole1))
    {
        camera1->projectPoint(-Transform12.trans(), epipole1);
        epipoleInverted1 = true;
    }
    else epipoleInverted1 = false;
    
    if (not camera2->projectPoint(Transform12.transInv(), epipole2))
    {
        camera2->projectPoint(-Transform12.transInv(), epipole2);
        epipoleInverted2 = false;
    }
    else epipoleInverted2 = false;
    
    epipolePx1 = round(epipole1);
    epipolePx2 = round(epipole2);
    
}

CurveRasterizer<int, Polynomial2> EnhancedStereo::getCurveRasteriser1(int idx)
{
    Vector2i pt = pointPxVec1[idx];
    CurveRasterizer<int, Polynomial2> raster(pt, epipolePx1, epipolar.getFirst(reconstVec[idx]));
    if (epipoleInverted1) raster.setStep(-1);
    return raster;
}

CurveRasterizer<int, Polynomial2> EnhancedStereo::getCurveRasteriser2(int idx)
{
    Vector2i pinfPx = pinfPxVec[idx];
    CurveRasterizer<int, Polynomial2> raster(pinfPx, epipolePx2, epipolar.getSecond(reconstVec[idx]));
    if (epipoleInverted2) raster.setStep(-1);
    return raster;
}

//TODO reconstruct the depth points, not everything
void EnhancedStereo::computeReconstructed()
{
    pointVec1.reserve(params.yMax*params.xMax);
    pointPxVec1.reserve(params.yMax*params.xMax);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            pointVec1.emplace_back(params.u(x), params.v(y));
            pointPxVec1.emplace_back(params.u(x), params.v(y));
        }
    }
    camera1->reconstructPointCloud(pointVec1, reconstVec, maskVec);
}

void EnhancedStereo::computeRotated()
{
    Transform12.inverseRotate(reconstVec, reconstRotVec);
}

void EnhancedStereo::computePinf()
{
    camera2->projectPointCloud(reconstRotVec, pinfVec);
    pinfPxVec.resize(pinfVec.size());
    for (int i = 0; i < pinfVec.size(); i++)
    {
        if (not maskVec[i]) continue;
        pinfPxVec[i] = round(pinfVec[i]);
    }
}

//TODO move elsewhere
void EnhancedStereo::traceEpipolarLine(int x, int y, Mat8u & out, CameraIdx camIdx)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::traceEpipolarLine" << endl;
    int idx = getLinearIndex(x, y);
    if (not maskVec[idx]) return;
    
    CurveRasterizer<int, Polynomial2> * raster = NULL;
    int count;
    if (camIdx == CAMERA_1)
    { 
        count = (pinfPxVec[idx] - epipolePx1).norm();
        raster = new CurveRasterizer<int, Polynomial2>(getCurveRasteriser1(idx));
    }
    else if (camIdx == CAMERA_2)
    {
        count = (pinfPxVec[idx] - epipolePx2).norm();
        raster = new CurveRasterizer<int, Polynomial2>(getCurveRasteriser2(idx));
    }
    else 
    {
        return;
    }
    for (int i = 0; i < count; i++)
    {
        out(raster->v, raster->u) = 0;
        out(raster->v+1, raster->u) = 0;
        out(raster->v, raster->u+1) = 0;
        out(raster->v+1, raster->u+1) = 0;
        raster->step();
    }
    delete raster;
}

void EnhancedStereo::createBuffer()
{
    if (params.verbosity > 1) cout << "EnhancedStereo::createBuffer" << endl;
    int bufferWidth = params.xMax*params.dispMax;
    if (errorBuffer.cols != bufferWidth or errorBuffer.rows != params.yMax)
    {
        errorBuffer = Mat8u(Size(bufferWidth, params.yMax));
    }
    if (tableauLeft.cols != bufferWidth or tableauLeft.rows != params.yMax)
    {
        tableauLeft = Mat32s(Size(bufferWidth, params.yMax));
    }
    if (tableauRight.cols != bufferWidth or tableauRight.rows != params.yMax)
    {
        tableauRight = Mat32s(Size(bufferWidth, params.yMax));
    }
    if (tableauTop.cols != bufferWidth or tableauTop.rows != params.yMax)
    {
        tableauTop = Mat32s(Size(bufferWidth, params.yMax));
    }
    if (tableauBottom.cols != bufferWidth or tableauBottom.rows != params.yMax)
    {
        tableauBottom = Mat32s(Size(bufferWidth, params.yMax));
    }
    if (smallDisparity.cols != params.xMax or smallDisparity.rows != params.yMax)
    {
        smallDisparity = Mat8u(Size(params.xMax, params.yMax));

    }
    if (params.verbosity > 2) 
    {
        cout << "    small disparity size: " << smallDisparity.size() << endl;
    }
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depth)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    depth = DepthMap(camera1, params.yMax, params.xMax,
            params.u0, params.v0, params.scale);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            computeDepth(x, y, depth.at(x, y), depth.sigma(x, y));
        }
    }
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, const Mat8u & img2, Mat32f & depthMat)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    computeDepth(depthMat);
}

//TODO change u,v to x,y
//TODO make step change scaleable and define in in a loop
void EnhancedStereo::computeCurveCost(const Mat8u & img1, const Mat8u & img2)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeCurveCost" << endl;
    
    int HALF_LENGTH = max(params.scale - 1, 1);
    int LENGTH = HALF_LENGTH * 2 + 1;
    
    // compute the weights for matching cost
    vector<int> kernelVec, waveVec;
    const int NORMALIZER = initKernel(kernelVec, LENGTH);
    const int WAVE_NORM = initWave(waveVec, LENGTH);
    EpipolarDescriptor epipolarDescriptor(LENGTH, WAVE_NORM, waveVec.data(), {1, 2, 3, 5, 8, 13, 21});
    vector<int> discHist(50, 0);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            if (params.verbosity > 3) 
            {
                cout << "    x: " << x << " y: " << y << "  idx: " << idx; 
                cout << "  mask: " << maskVec[idx] <<  endl;
            }
            if (not maskVec[idx])
            {
                uint8_t * outPtr = errorBuffer.row(y).data + x*params.dispMax;            
                *outPtr = 0;
                fill(outPtr + 1, outPtr + params.dispMax, 255);
                continue;
            }
            // compute the local image descriptor -- a piece of the epipolar curve on the first image
            vector<uint8_t> descriptor;
            CurveRasterizer<int, Polynomial2> descRaster = getCurveRasteriser1(idx);
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            discHist[step]++;
            int nSteps = (params.dispMax  + step - 1 ) / step; 
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
               
            //sample the curve 
            raster.setStep(step); 
            raster.steps(-HALF_LENGTH);           
            vector<uint8_t> sampleVec(nSteps + LENGTH - 1, 0);
            for (int i = 0; i  < nSteps + LENGTH - 1; i++, raster.step())
            {
                if (raster.v < 0 or raster.v >= img2.rows 
                    or raster.u < 0 or raster.u >= img2.cols) sampleVec[i] = 0;
                else sampleVec[i] = img2(raster.v, raster.u);
            }
            
            //compute the bias;
            int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
            
            // fill up the cost buffer
            uint8_t * outPtr = errorBuffer.row(y).data + x*params.dispMax;
            for (int d = 0; d < nSteps; d++, outPtr += step)
            {
                int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
                int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
                int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
                                descriptor.begin(), sampleVec.begin() + d, bias);
                *outPtr = acc / NORMALIZER;
            }
            if (step > 1) fillGaps(errorBuffer.row(y).data + x*params.dispMax, step);
//            }
//            else
//            {
//                //sample img2 along the epipolar curve
//                
//                raster.steps(-HALF_LENGTH);
//                vector<uint8_t> sampleVec(params.dispMax + LENGTH - 1, 0);
//                for (int i = 0; i  < params.dispMax + LENGTH - 1; i++, raster.step())
//                {
//                    if (raster.v < 0 or raster.v >= img2.rows 
//                        or raster.u < 0 or raster.u >= img2.cols) sampleVec[i] = 0;
//                    else sampleVec[i] = img2(raster.v, raster.u);
//                }
//                
//                //compute the bias;
//                int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
//                // fill up the cost buffer
//                uint8_t * outPtr = errorBuffer.row(y).data + x*params.dispMax;
//                for (int d = 0; d < params.dispMax; d++, outPtr++)
//                {
//                    int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                    int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
//                    int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
//                                    descriptor.begin(), sampleVec.begin() + d, bias);
//                    *outPtr = acc / NORMALIZER;
//                }
//            }
            
        }
        
    }
}

void EnhancedStereo::fillGaps(uint8_t * const data, const int step)
{
    assert(step > 0);
    //linear interpolation for all intermediate points
    int base;
    switch (step)
    {
    case 2:
        for (base = 2; base < params.dispMax; base += 2)
        {
            data[base - 1] = (data[base - 2] + data[base]) >> 1;
        }
        break;
    case 3:
        for (base = 3; base < params.dispMax; base += 3)
        {
            const uint8_t & val1 = data[base - 3];
            const uint8_t & val2 = data[base];
            data[base - 2] = (val1 << 1  + val2) / 3;
            data[base - 1] = (val1 + val2 << 1) / 3;
        }
        break;
    default:
        for (base = step; base < params.dispMax; base += step)
        {
            const uint8_t & val1 = data[base - step];
            const uint8_t & val2 = data[base];
            for (int i = step - 1; i > 0; i--)
            {
                data[base - i] = (val1 * i + val2 * (step - i)) / step;
            }
        }
        break;
    }
    //for the rest just constant extrapolation
    base -= step;
    const uint8_t & val = data[base];
    for (int i = base + 1 ; i < params.dispMax; i++)
    {
        data[i] = val;
    }
}

void EnhancedStereo::computeDynamicStep(const int32_t* inCost, const uint8_t * error, int32_t * outCost)
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
    for (int y = 0; y < params.yMax; y++)
    {
        int32_t * tableauRow = (int32_t *)(tableauLeft.row(y).data);
        uint8_t * errorRow = errorBuffer.row(y).data;
        // init the first row
        copy(errorRow, errorRow + params.dispMax, tableauRow);
        // fill up the tableau
        for (int x = 1; x < params.xMax; x++)
        {
            computeDynamicStep(tableauRow + (x - 1)*params.dispMax,
                    errorRow + x*params.dispMax, tableauRow + x*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    right" << endl;  
    // right tableau init
    for (int y = 0; y < params.yMax; y++)
    {
        int32_t * tableauRow = (int32_t *)(tableauRight.row(y).data);
        uint8_t * errorRow = errorBuffer.row(y).data;
        int base = (params.xMax - 1) * params.dispMax;
        copy(errorRow + base, errorRow + base + params.dispMax, tableauRow + base);
        
        for (int x = params.xMax - 2; x >= 0; x--)
        {
            computeDynamicStep(tableauRow + (x + 1)*params.dispMax, 
                    errorRow + x*params.dispMax, tableauRow + x*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    top" << endl;
    // top-down tableau init
    for (int x = 0; x < params.xMax; x++)
    {
        auto tableauCol = tableauTop(Rect(x*params.dispMax, 0, params.dispMax, params.yMax));
        auto errorCol = errorBuffer(Rect(x*params.dispMax, 0, params.dispMax, params.yMax));
        copy(errorCol.data, errorCol.data + params.dispMax, (int*)(tableauCol.data));
        for (int y = 1; y < params.yMax; y++)
        {
            computeDynamicStep((int32_t*)(tableauCol.row(y-1).data), 
                    errorCol.row(y).data,
                    (int32_t*)(tableauCol.row(y).data));
        }
    }
    if (params.verbosity > 1) cout << "    bottom" << endl;
    // bottom-up tableau init
    for (int x = 0; x < params.xMax; x++)
    {
        auto tableauCol = tableauBottom(Rect(x*params.dispMax, 0, params.dispMax, params.yMax));
        auto errorCol = errorBuffer(Rect(x*params.dispMax, 0, params.dispMax, params.yMax));
        int vLast = params.yMax - 1;
        copy(errorCol.row(vLast).data, 
                errorCol.row(vLast).data + params.dispMax, 
                (int32_t*)(tableauCol.row(vLast).data));
        for (int y = params.yMax - 2; y >= 0; y--)
        {
            computeDynamicStep((int32_t*)(tableauCol.row(y+1).data), 
                    errorCol.row(y).data,
                    (int32_t*)(tableauCol.row(y).data));
        }
    }
    
}

void EnhancedStereo::reconstructDisparity()
{
    if (params.verbosity > 0) cout << "EnhancedStereo::reconstructDisparity" << endl;
    Mat16s errFinalMat(smallDisparity.size());
    for (int y = 0; y < params.yMax; y++)
    {
        int32_t* dynRow1 = (int32_t*)(tableauLeft.row(y).data);
        int32_t* dynRow2 = (int32_t*)(tableauRight.row(y).data);
        int32_t* dynRow3 = (int32_t*)(tableauTop.row(y).data);
        int32_t* dynRow4 = (int32_t*)(tableauBottom.row(y).data);
        uint8_t* errRow = (errorBuffer.row(y).data);
        for (int x = 0; x < params.xMax; x++)
        {
            int bestCost = 1000000;
            uint8_t & bestDisp = smallDisparity(y, x);
            int16_t & bestErr = errFinalMat(y, x);
            bestDisp = 0;
            for (int d = 0; d < params.dispMax; d++)
            {
                int base = x * params.dispMax;
                int acc = dynRow1[base + d] + dynRow2[base + d] 
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
    if (not camera1->reconstructPoint(Vector2d(x1, y1), v1) or 
        not camera2->reconstructPoint(Vector2d(x2, y2), v2) )
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

void EnhancedStereo::computeDepth(Mat32f & distance)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeDepth(Mat32f &)" << endl;
    distance.create(params.yMax, params.xMax);
    for (int v = 0; v < params.yMax; v++)
    {
        for (int u = 0; u < params.xMax; u++)
        {
            distance(v, u) = computeDepth(u, v);
        }
    }
}

bool EnhancedStereo::computeDepth(int x, int y, double & dist, double & sigma)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::computeDepth(int *)" << endl;
    
    int idx = getLinearIndex(x, y);
    if (not maskVec[idx])
    { 
        dist = 0;
        sigma = 0;
        return false;
    
    }
    int disparity = smallDisparity(y, x);
    if (disparity <= 0) 
    {
        dist = params.maxDepth;
        sigma = params.maxDepth;
        return true;
    }
    
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    // to compute point on the second image
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
    raster.steps(disparity - 1);
    
    
    Vector3d X;
    if (triangulate(pt1[0], pt1[1], raster.u, raster.v, X))
    {
        sigma = X.norm();
        raster.step();
        if (triangulate(pt1[0], pt1[1], raster.u, raster.v, X))
        {
            dist = X.norm();
            sigma = abs(sigma - dist) / 2;
        }
        else
        {
            sigma = 0;
            return false;
        }
    }
    else 
    {
        dist = 0;
        sigma = 0;
        return false;
    }
}

double EnhancedStereo::computeDepth(int u, int v)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::computeDepth(int *)" << endl;
    
    int idx = getLinearIndex(u, v);
    if (not maskVec[idx]) return 0;
    int disparity = smallDisparity(v, u);
    if (disparity <= 0) 
    {
        return params.maxDepth;
    }
    
    // to compute point on the second image
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
    raster.steps(disparity);
    
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    Vector3d X;
    if (triangulate(pt1[0], pt1[1], raster.u, raster.v, X))
    {
        return X.norm();
    }
    else 
    {
        return 0;
    }
}

//TODO remove this function
void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        Mat32f & distanceMat, const Vector3dVec & polygonVec)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::generatePlane" << endl;
    distanceMat.create(params.yMax, params.xMax);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            distanceMat(y, x) = 0;
            Vector3d vec; // the direction vector
            if (not camera1->reconstructPoint(Vector2d(params.u(x), params.v(y)), vec)) continue;
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
            distanceMat(y, x) = vec.norm();
        }
    }
}

void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        DepthMap & depth, const Vector3dVec & polygonVec)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::generatePlane" << endl;
    depth = DepthMap(camera1, params.xMax, params.yMax, params.u0, params.v0, params.scale);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int v = 0; v < params.yMax; v++)
    {
        for (int u = 0; u < params.xMax; u++)
        {
            depth.at(u, v) = 0;
            Vector3d vec; // the direction vector
            if (not camera1->reconstructPoint(Vector2d(params.u(u), params.v(v)), vec)) continue;
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
            depth.at(u, v) = vec.norm();
        }
    }
}

void EnhancedStereo::upsampleDisparity(const Mat8u & img1, Mat8u & disparityMat)
{
    cout << smallDisparity.size() << endl;
    smallDisparity.copyTo(disparityMat);
//    resize(smallDisparity, disparity, Size(0, 0), params.scale, params.scale, 0);
}

