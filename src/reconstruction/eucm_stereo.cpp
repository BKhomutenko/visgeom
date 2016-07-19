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
    if (not cam1.projectPoint(Transform12.trans(), epipole1))
    {
        cam1.projectPoint(-Transform12.trans(), epipole1);
        epipoleInverted1 = true;
    }
    else epipoleInverted1 = false;
    
    if (not cam2.projectPoint(Transform12.transInv(), epipole2))
    {
        cam2.projectPoint(-Transform12.transInv(), epipole2);
        epipoleInverted2 = false;
    }
    else epipoleInverted2 = false;
    
    for (int i = 0; i < 2; i++)
    {
        epipolePx1[i] = round(epipole1[i]);
        epipolePx2[i] = round(epipole2[i]);
    }
    
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
    pointVec1.reserve(params.dispHeight*params.dispWidth);
    pointPxVec1.reserve(params.dispHeight*params.dispWidth);
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            pointVec1.emplace_back(params.uImg(u), params.vImg(v));
            pointPxVec1.emplace_back(params.uImg(u), params.vImg(v));
        }
    }
    cam1.reconstructPointCloud(pointVec1, reconstVec, maskVec);
}

// computes epipolarDirectionVec -- computed by shifting the reconstructed points in the direction 
// of motion infinitesimally and projecting them back
// TODO replace by epipolar rasterizer
void EnhancedStereo::computeEpipolarDirections()
{
    epipolarDirectionVec.resize(pointVec1.size());
    Vector3d t = Transform12.trans();
    t = t.normalized() * 0.001;
    for (int i = 0; i < pointVec1.size(); i++)
    {
        if (not maskVec[i]) continue;
        Vector3d Xt = reconstVec[i] - t;
        Vector2d pt;
        cam1.projectPoint(Xt, pt);
        pt -= pointVec1[i];
        epipolarDirectionVec[i] = pt.normalized();
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
        if (not maskVec[i]) continue;
        pinfPxVec[i][0] = round(pinfVec[i][0]);
        pinfPxVec[i][1] = round(pinfVec[i][1]);
    }
}

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
        out(raster->y, raster->x) = 0;
        out(raster->y+1, raster->x) = 0;
        out(raster->y, raster->x+1) = 0;
        out(raster->y+1, raster->x+1) = 0;
        raster->step();
    }
    delete raster;
}

void EnhancedStereo::createBuffer()
{
    if (params.verbosity > 1) cout << "EnhancedStereo::createBuffer" << endl;
    int bufferWidth = params.dispWidth*params.dispMax;
    if (errorBuffer.cols != bufferWidth or errorBuffer.rows != params.dispHeight)
    {
        errorBuffer = Mat8u(Size(bufferWidth, params.dispHeight));
    }
    if (tableauLeft.cols != bufferWidth or tableauLeft.rows != params.dispHeight)
    {
        tableauLeft = Mat32s(Size(bufferWidth, params.dispHeight));
    }
    if (tableauRight.cols != bufferWidth or tableauRight.rows != params.dispHeight)
    {
        tableauRight = Mat32s(Size(bufferWidth, params.dispHeight));
    }
    if (tableauTop.cols != bufferWidth or tableauTop.rows != params.dispHeight)
    {
        tableauTop = Mat32s(Size(bufferWidth, params.dispHeight));
    }
    if (tableauBottom.cols != bufferWidth or tableauBottom.rows != params.dispHeight)
    {
        tableauBottom = Mat32s(Size(bufferWidth, params.dispHeight));
    }
    if (smallDisparity.cols != params.dispWidth or smallDisparity.rows != params.dispHeight)
    {
        smallDisparity = Mat8u(Size(params.dispWidth, params.dispHeight));

    }
    if (params.verbosity > 2) 
    {
        cout << "    small disparity size: " << smallDisparity.size() << endl;
    }
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, 
        const Mat8u & img2,
        Mat8u & disparityMat)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    smallDisparity.copyTo(disparityMat);
}

void EnhancedStereo::comuteStereo(const Mat8u & img1, 
        const Mat8u & img2,
        DepthMap & depth)
{
    computeCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    //TODO fuse the following code with computeDistance
    depth = DepthMap(&cam1, params.dispHeight, params.dispWidth,
            params.u0, params.v0, params.scale);
    for (int y = 0; y < params.dispHeight; y++)
    {
        for (int x = 0; x < params.dispWidth; x++)
        {
            depth.at(x, y) = computeDistance(x, y);
        }
    }
}


//TODO change u,v to x,y
//TODO possibly change the step length along the curve to reduce the computation cost
void EnhancedStereo::computeCurveCost(const Mat8u & img1, const Mat8u & img2)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::computeCurveCost" << endl;
    
    int HALF_LENGTH = max(params.scale - 1, 1);
    int LENGTH = HALF_LENGTH * 2 + 1;
    
    // compute the weights for matching cost
    vector<int> kernelVec(LENGTH);
    vector<int> waveVec(LENGTH);
    int WAVE_NORM, NORMALIZER;
    switch (LENGTH)
    {
    case 3:
        copy(KERNEL_3.begin(), KERNEL_3.end(), kernelVec.begin());
        copy(WAVE_3.begin(), WAVE_3.end(), waveVec.begin());
        NORMALIZER = NORMALIZER_3;
        WAVE_NORM = WAVE_NORM_3;
        break;
    default:
        LENGTH = 5;
        HALF_LENGTH = 2;
    case 5:
        copy(KERNEL_5.begin(), KERNEL_5.end(), kernelVec.begin());
        copy(WAVE_5.begin(), WAVE_5.end(), waveVec.begin());
        NORMALIZER = NORMALIZER_5;
        WAVE_NORM = WAVE_NORM_5;
        break;
    case 7:
        copy(KERNEL_7.begin(), KERNEL_7.end(), kernelVec.begin());
        copy(WAVE_7.begin(), WAVE_7.end(), waveVec.begin());
        NORMALIZER = NORMALIZER_7;
        WAVE_NORM = WAVE_NORM_7;
        break;
    case 9:
        copy(KERNEL_9.begin(), KERNEL_9.end(), kernelVec.begin());
        copy(WAVE_9.begin(), WAVE_9.end(), waveVec.begin());
        NORMALIZER = NORMALIZER_9;
        WAVE_NORM = WAVE_NORM_9;
        break;
    }
    WAVE_NORM *= 2;
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            int idx = getLinearIndex(u, v);
            if (params.verbosity > 3) 
            {
                cout << "    x: " << u << " y: " << v << "  idx: " << idx; 
                cout << "  mask: " << maskVec[idx] <<  endl;
            }
            if (not maskVec[idx])
            {
                uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;            
                *outPtr = 0;
                fill(outPtr + 1, outPtr + params.dispMax, 255);
                continue;
            }
            // compute the local image descriptor -- a piece of the epipolar curve on the first image
            CurveRasterizer<int, Polynomial2> descRaster = getCurveRasteriser1(idx);
            descRaster.setStep(-1);
            descRaster.steps(-HALF_LENGTH);
            vector<uint8_t> descriptor(LENGTH);
            for (int i = 0; i < LENGTH; i++, descRaster.step())
            {
                descriptor[i] = img1(descRaster.y, descRaster.x);
            }
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
            int descAcc = 0;
            for (int i = 0, k = 1; i < LENGTH; i++)
            {
                descAcc += waveVec[i] *  descriptor[i];
            }
            if (abs(descAcc) < WAVE_NORM) //TODO check the constant
            {
                
                raster.setStep(2);  
                
//                recompute the descriptor
                descRaster = getCurveRasteriser1(idx);
                descRaster.setStep(-2);
                descRaster.steps(-HALF_LENGTH);
                for (int i = 0; i < LENGTH; i++, descRaster.step())
                {
                    descriptor[i] = img1(descRaster.y, descRaster.x);
                }
            
                raster.steps(-HALF_LENGTH);
                
                vector<uint8_t> sampleVec(params.dispMax/2 + LENGTH - 1, 0);
                for (int i = 0; i  < params.dispMax/2 + LENGTH - 1; i++, raster.step())
                {
                    if (raster.y < 0 or raster.y >= img2.rows 
                        or raster.x < 0 or raster.x >= img2.cols) sampleVec[i] = 0;
                    else sampleVec[i] = img2(raster.y, raster.x);
                }
                
                //compute the bias;
                int sum1 = accumulate(descriptor.begin(), descriptor.end(), 0);
                
                // fill up the cost buffer
                uint8_t * outPtr = errorBuffer.row(v).data + u*params.dispMax;
                for (int d = 0; d < params.dispMax/2; d++, outPtr+=2)
                {
                    int acc = 0;
                    int sum2 = accumulate(sampleVec.begin() + d, sampleVec.begin() + d + LENGTH, 0);
                    int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
                    for (int i = 0; i < LENGTH; i++)
                    {
                        acc += abs(descriptor[i] - sampleVec[d + i] + bias) * kernelVec[i];
                    }
                    *outPtr = acc / NORMALIZER;
                    if (d)
                    {
                        *(outPtr - 1) = (*(outPtr - 2) >> 1) + (*outPtr >> 1);
                    }
                }
                *(outPtr - 1) = *(outPtr - 2);
                
            }
            else
            {
                //sample img2 along the epipolar curve
                
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
                        acc += abs(descriptor[i] - sampleVec[d + i] + bias) * kernelVec[i];
                    }
                    
                    *outPtr = acc / NORMALIZER;
                }
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
    
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            int idx = getLinearIndex(u, v);
            if (not maskVec[idx]) continue;
            //FIXME this makes the algorithm dependent on the motion direction
            
            
//            Point pPxinf = pinfPxVec[idx];
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
            
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
            int u1 = params.uImg(u) - hblock;
            int v1 = params.vImg(v) - hblock;
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
    for (int v = 0; v < params.dispHeight; v++)
    {
        int32_t * tableauRow = (int32_t *)(tableauLeft.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        // init the first row
        copy(errorRow, errorRow + params.dispMax, tableauRow);
        // fill up the tableau
        for (int u = 1; u < params.dispWidth; u++)
        {
            computeDynamicStep(tableauRow + (u - 1)*params.dispMax,
                    errorRow + u*params.dispMax, tableauRow + u*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    right" << endl;  
    // right tableau init
    for (int v = 0; v < params.dispHeight; v++)
    {
        int32_t * tableauRow = (int32_t *)(tableauRight.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        int base = (params.dispWidth - 1) * params.dispMax;
        copy(errorRow + base, errorRow + base + params.dispMax, tableauRow + base);
        
        for (int u = params.dispWidth - 2; u >= 0; u--)
        {
            computeDynamicStep(tableauRow + (u + 1)*params.dispMax, 
                    errorRow + u*params.dispMax, tableauRow + u*params.dispMax);
        }
    }
    if (params.verbosity > 1) cout << "    top" << endl;
    // top-down tableau init
    for (int u = 0; u < params.dispWidth; u++)
    {
        auto tableauCol = tableauTop(Rect(u*params.dispMax, 0, params.dispMax, params.dispHeight));
        auto errorCol = errorBuffer(Rect(u*params.dispMax, 0, params.dispMax, params.dispHeight));
        copy(errorCol.data, errorCol.data + params.dispMax, (int*)(tableauCol.data));
        for (int v = 1; v < params.dispHeight; v++)
        {
            computeDynamicStep((int32_t*)(tableauCol.row(v-1).data), 
                    errorCol.row(v).data,
                    (int32_t*)(tableauCol.row(v).data));
        }
    }
    if (params.verbosity > 1) cout << "    bottom" << endl;
    // bottom-up tableau init
    for (int u = 0; u < params.dispWidth; u++)
    {
        auto tableauCol = tableauBottom(Rect(u*params.dispMax, 0, params.dispMax, params.dispHeight));
        auto errorCol = errorBuffer(Rect(u*params.dispMax, 0, params.dispMax, params.dispHeight));
        int vLast = params.dispHeight - 1;
        copy(errorCol.row(vLast).data, 
                errorCol.row(vLast).data + params.dispMax, 
                (int32_t*)(tableauCol.row(vLast).data));
        for (int v = params.dispHeight - 2; v >= 0; v--)
        {
            computeDynamicStep((int32_t*)(tableauCol.row(v+1).data), 
                    errorCol.row(v).data,
                    (int32_t*)(tableauCol.row(v).data));
        }
    }
    
}

void EnhancedStereo::reconstructDisparity()
{
    if (params.verbosity > 0) cout << "EnhancedStereo::reconstructDisparity" << endl;
    Mat16s errFinalMat(smallDisparity.size());
    for (int v = 0; v < params.dispHeight; v++)
    {
        int32_t* dynRow1 = (int32_t*)(tableauLeft.row(v).data);
        int32_t* dynRow2 = (int32_t*)(tableauRight.row(v).data);
        int32_t* dynRow3 = (int32_t*)(tableauTop.row(v).data);
        int32_t* dynRow4 = (int32_t*)(tableauBottom.row(v).data);
        uint8_t* errRow = (errorBuffer.row(v).data);
        for (int u = 0; u < params.dispWidth; u++)
        {
            int bestCost = 1000000;
            uint8_t & bestDisp = smallDisparity(v, u);
            int16_t & bestErr = errFinalMat(v, u);
            bestDisp = 0;
            for (int d = 0; d < params.dispMax; d++)
            {
                int base = u * params.dispMax;
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
    distance.create(params.dispHeight, params.dispWidth);
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            distance(v, u) = computeDistance(u, v);
        }
    }
}

bool EnhancedStereo::computeDistance(int x, int y, double & dist, double & sigma)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::computeDistance(int *)" << endl;
    
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
        dist = params.maxDistance;
        sigma = params.maxDistance;
        return true;
    }
    
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    // to compute point on the second image
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
    raster.steps(disparity - 1);
    
    
    Vector3d X;
    if (triangulate(pt1[0], pt1[1], raster.x, raster.y, X))
    {
        sigma = X.norm();
        raster.step();
        if (triangulate(pt1[0], pt1[1], raster.x, raster.y, X))
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

double EnhancedStereo::computeDistance(int u, int v)
{
    if (params.verbosity > 3) cout << "EnhancedStereo::computeDistance(int *)" << endl;
    
    int idx = getLinearIndex(u, v);
    if (not maskVec[idx]) return 0;
    int disparity = smallDisparity(v, u);
    if (disparity <= 0) 
    {
        return params.maxDistance;
    }
    
    // to compute point on the second image
    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
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

//TODO remove this function
void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        Mat32f & distanceMat, const Vector3dVec & polygonVec)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::generatePlane" << endl;
    distanceMat.create(params.dispHeight, params.dispWidth);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            distanceMat(v, u) = 0;
            Vector3d vec; // the direction vector
            if (not cam1.reconstructPoint(Vector2d(params.uImg(u), params.vImg(v)), vec)) continue;
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
            distanceMat(v, u) = vec.norm();
        }
    }
}

void EnhancedStereo::generatePlane(Transformation<double> TcameraPlane,
        DepthMap & depthMap, const Vector3dVec & polygonVec)
{
    if (params.verbosity > 0) cout << "EnhancedStereo::generatePlane" << endl;
    depthMap = DepthMap(&cam1, params.dispWidth, params.dispHeight, params.u0, params.v0, params.scale);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int v = 0; v < params.dispHeight; v++)
    {
        for (int u = 0; u < params.dispWidth; u++)
        {
            depthMap.at(u, v) = 0;
            Vector3d vec; // the direction vector
            if (not cam1.reconstructPoint(Vector2d(params.uImg(u), params.vImg(v)), vec)) continue;
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
            depthMap.at(u, v) = vec.norm();
        }
    }
}

void EnhancedStereo::upsampleDisparity(const Mat8u & img1, Mat8u & disparityMat)
{
    cout << smallDisparity.size() << endl;
    smallDisparity.copyTo(disparityMat);
//    resize(smallDisparity, disparity, Size(0, 0), params.scale, params.scale, 0);
}

