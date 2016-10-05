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
#include "utils/filter.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

CurveRasterizer<int, Polynomial2> EnhancedSGM::getCurveRasteriser1(int idx)
{
    Vector2i pt = pointPxVec1[idx];
    CurveRasterizer<int, Polynomial2> raster(pt, epipoles.getFirstPx(), epipolar.getFirst(reconstVec[idx]));
    if (epipoles.firstIsInverted()) raster.setStep(-1);
    return raster;
}

CurveRasterizer<int, Polynomial2> EnhancedSGM::getCurveRasteriser2(int idx)
{
    Vector2i pinfPx = pinfPxVec[idx];
    CurveRasterizer<int, Polynomial2> raster(pinfPx, epipoles.getSecondPx(),
             epipolar.getSecond(reconstVec[idx]));
    if (epipoles.secondIsInverted()) raster.setStep(-1);
    return raster;
}

//TODO reconstruct the depth points, not everything
void EnhancedSGM::computeReconstructed()
{
    pointVec1.reserve(params.yMax*params.xMax);
    pointPxVec1.reserve(params.yMax*params.xMax);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            pointVec1.emplace_back(params.uConv(x), params.vConv(y));
            pointPxVec1.emplace_back(params.uConv(x), params.vConv(y));
        }
    }
    camera1->reconstructPointCloud(pointVec1, reconstVec, maskVec);
}

void EnhancedSGM::computeRotated()
{
    Transform12.inverseRotate(reconstVec, reconstRotVec);
}

//FIXME maskVec must be recomputed to discard not projected pInf
void EnhancedSGM::computePinf()
{
    camera2->projectPointCloud(reconstRotVec, pinfVec);
    pinfPxVec.resize(pinfVec.size());
    for (int i = 0; i < pinfVec.size(); i++)
    {
        if (not maskVec[i]) continue;
        pinfPxVec[i] = round(pinfVec[i]);
    }
}

void EnhancedSGM::computeUVCache()
{
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
            raster.steps(-DISPARITY_MARGIN);
            const int uvCacheStep = params.dispMax + 2 * DISPARITY_MARGIN;
            int32_t * uPtr = (int32_t *)uCache.row(y).data + x*uvCacheStep;
            int32_t * vPtr = (int32_t *)vCache.row(y).data + x*uvCacheStep;
            for (int i = 0; i  < uvCacheStep; i++, raster.step(), uPtr++, vPtr++)
            {
                if (raster.v < 0 or raster.v >= params.vMax 
                    or raster.u < 0 or raster.u >= params.uMax)
                {
                    // coordinate is out of the image
                    *uPtr = -1;
                    *vPtr = -1;
                }
                else
                {
                    // coordinate is within the image
                    *uPtr = raster.u;
                    *vPtr = raster.v;
                }
            }
        }
    }
}

void EnhancedSGM::createBuffer()
{
    if (params.verbosity > 1) cout << "EnhancedSGM::createBuffer" << endl;
    assert(params.hypMax > 0);
    int bufferWidth = params.xMax*params.dispMax;
    errorBuffer.create(params.yMax, bufferWidth);
    tableauLeft.create(params.yMax, bufferWidth);
    tableauRight.create(params.yMax, bufferWidth);
    tableauTop.create(params.yMax, bufferWidth);
    tableauBottom.create(params.yMax, bufferWidth);
    smallDisparity.create(params.yMax, params.xMax * params.hypMax);
    finalErrorMat.create(params.yMax, params.xMax * params.hypMax);
    if (params.imageBasedCost) costBuffer.create(params.yMax, params.xMax);
    if (params.salientPoints) salientBuffer.create(params.yMax, params.xMax);
    uCache.create(params.yMax, params.xMax * (params.dispMax + 2*DISPARITY_MARGIN));
    vCache.create(params.yMax, params.xMax * (params.dispMax + 2*DISPARITY_MARGIN));
    if (params.verbosity > 2) 
    {
        cout << "    small disparity size: " << smallDisparity.size() << endl;
    }
}

void EnhancedSGM::computeStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depth)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    depth = DepthMap(camera1, params, params.hypMax);
    assert((ScaleParameters)depth == (ScaleParameters)params);
    for (int h = 0; h < params.hypMax; h++)
    {
        for (int y = 0; y < params.yMax; y++)
        {
            for (int x = 0; x < params.xMax; x++)
            {
                if (not params.salientPoints or salientBuffer(y, x)) 
                {
                    //FIXME temporary 
//                    depth.at(x, y, h) =  smallDisparity(y, x*params.hypMax + h);
                    computeDepth(depth.at(x, y, h), depth.sigma(x, y, h), x, y, h);
                    depth.cost(x, y, h) = errorBuffer(y, x*params.dispMax + h);
                }
                else 
                {
                    depth.at(x, y, h) = OUT_OF_RANGE;
                    depth.sigma(x, y, h) = OUT_OF_RANGE;
                    depth.cost(x, y, h) = 100;
                }
            }
        }
    }
}

//TODO remove, depricated function
void EnhancedSGM::computeStereo(const Mat8u & img1, const Mat8u & img2, Mat32f & depthMat)
{
    computeCurveCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    computeDepth(depthMat);
}

void EnhancedSGM::computeCurveCost(const Mat8u & img1, const Mat8u & img2)
{
    if (params.verbosity > 0) cout << "EnhancedSGM::computeCurveCost" << endl;
    
    
    const int HALF_LENGTH = getHalfLength();
    const int LENGTH = HALF_LENGTH * 2 + 1;
    cout << "SGM LENGTH: " << LENGTH << endl;
    // compute the weights for matching cost
    vector<int> kernelVec, waveVec;
    const int NORMALIZER = initKernel(kernelVec, LENGTH);
    const int WAVE_NORM = initWave(waveVec, LENGTH);
    EpipolarDescriptor epipolarDescriptor(LENGTH, 10, waveVec.data(), {1, 2, 3, 5});
    
    if (params.salientPoints) salientBuffer.setTo(0);
    
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            if (params.verbosity > 4) 
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
            // compute the local image descriptor,
            // a piece of the epipolar curve on the first image
            vector<uint8_t> descriptor;
            CurveRasterizer<int, Polynomial2> descRaster = getCurveRasteriser1(idx);
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            if (step < 1) 
            {
                //TODO make a function
                uint8_t * outPtr = errorBuffer.row(y).data + x*params.dispMax;            
                *outPtr = 0;
                fill(outPtr + 1, outPtr + params.dispMax, 255);
                continue;
            }
            
            if (params.imageBasedCost) 
            {
                switch (step)
                {
                case 1:
                    costBuffer(y, x) = params.lambdaJump;
                    break;
                case 2:
                    costBuffer(y, x) = params.lambdaJump * 3;
                    break;
                default:
                    costBuffer(y, x) = params.lambdaJump * 6;
                    break;
                }
            }
            
            //TODO revise the criterion (step == 1)
            if (params.salientPoints and step <= 2)
            {
                salientBuffer(y, x) = 1;
            }
            
            const int nSteps = ( params.dispMax  + step - 1 ) / step; 
               
            //sample the curve 
            vector<uint8_t> sampleVec(nSteps + LENGTH - 1, 0);
            
            if (params.useUVCache)
            {
                const int uvCacheStep = params.dispMax + 2 * DISPARITY_MARGIN;
                int32_t * uPtr = (int32_t *)uCache.row(y).data + x*uvCacheStep;
                int32_t * vPtr = (int32_t *)vCache.row(y).data + x*uvCacheStep;
                uPtr += DISPARITY_MARGIN - HALF_LENGTH * step;
                vPtr += DISPARITY_MARGIN - HALF_LENGTH * step;
                for (int i = 0; i  < nSteps + LENGTH - 1; i++, uPtr += step, vPtr += step)
                {
                    if (*uPtr < 0 or *vPtr < 0) sampleVec[i] = 0;
                    else sampleVec[i] = img2(*vPtr, *uPtr);
                }
            }
            else
            {
                CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
                raster.setStep(step); 
                raster.steps(-HALF_LENGTH);           
                
                for (int i = 0; i  < nSteps + LENGTH - 1; i++, raster.step())
                {
                    if (raster.v < 0 or raster.v >= img2.rows 
                        or raster.u < 0 or raster.u >= img2.cols) sampleVec[i] = 0;
                    else sampleVec[i] = img2(raster.v, raster.u);
                }
            }
            vector<int> costVec = compareDescriptor(descriptor, sampleVec);
//            //compute the bias;
//            int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
            
            // fill up the cost buffer
            uint8_t * outPtr = errorBuffer.row(y).data + x*params.dispMax;
            auto costIter = costVec.begin() + HALF_LENGTH;
            for (int d = 0; d < nSteps; d++, outPtr += step)
            {
//                int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
//                int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
//                                descriptor.begin(), sampleVec.begin() + d, bias);
//                *outPtr = acc / NORMALIZER;

                *outPtr = *costIter / LENGTH;
                ++costIter;
            }
            if (step > 1) fillGaps(errorBuffer.row(y).data + x*params.dispMax, step);
        }
    }
}

void EnhancedSGM::fillGaps(uint8_t * const data, const int step)
{
    assert(step > 0);
    //linear interpolation for all intermediate points
    int base;
    switch (step)
    {
    case 2:
        for (base = 2; base < params.dispMax; base += 2)
        {
            data[base - 1] = (data[base - 2] + data[base]) / 2;
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

void EnhancedSGM::computeDynamicStep(const int32_t* inCost, const uint8_t * error, int32_t * outCost)
{
    int bestCost = inCost[0];
    for (int i = 1; i < params.dispMax; i++)
    {
        bestCost = min(bestCost, inCost[i]);
    }
    int & val0 = outCost[0];
    val0 = inCost[0];
    val0 = min(val0, inCost[1] + params.lambdaStep);
    val0 = min(val0, bestCost + jumpCost);
    val0 += error[0];
    for (int i = 1; i < params.dispMax-1; i++)
    {
        int & val = outCost[i];
        val = inCost[i];
        val = min(val, inCost[i + 1] + params.lambdaStep);
        val = min(val, inCost[i - 1] + params.lambdaStep);
        val = min(val, bestCost + jumpCost);
        val += error[i];
    }
    int & vald = outCost[params.dispMax - 1];
    vald = inCost[params.dispMax - 1];
    vald = min(vald, inCost[params.dispMax - 2] + params.lambdaStep);
    vald = min(vald, bestCost + jumpCost);
    vald += error[params.dispMax - 1];
}

void EnhancedSGM::computeDynamicProgramming()
{
    if (params.verbosity > 0) cout << "EnhancedSGM::computeDynamicProgramming" << endl;
    if (params.verbosity > 1) cout << "    left" << endl;
    
    if (not params.imageBasedCost) jumpCost = params.lambdaJump;
    
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
            if (params.imageBasedCost) jumpCost = costBuffer(y, x);
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
            if (params.imageBasedCost) jumpCost = costBuffer(y, x);
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
            if (params.imageBasedCost) jumpCost = costBuffer(y, x);
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
            if (params.imageBasedCost) jumpCost = costBuffer(y, x);
            computeDynamicStep((int32_t*)(tableauCol.row(y+1).data), 
                    errorCol.row(y).data,
                    (int32_t*)(tableauCol.row(y).data));
        }
    }
    
}

//TODO make with local minima
void EnhancedSGM::reconstructDisparity()
{
    if (params.verbosity > 0) cout << "EnhancedSGM::reconstructDisparity" << endl;
    const int hypShift = params.xMax*params.yMax;
//    int sizeAcc = 0;
//    int sizeCount = 0;
    for (int y = 0; y < params.yMax; y++)
    {
        int32_t* dynRow1 = (int32_t*)(tableauLeft.row(y).data);
        int32_t* dynRow2 = (int32_t*)(tableauRight.row(y).data);
        int32_t* dynRow3 = (int32_t*)(tableauTop.row(y).data);
        int32_t* dynRow4 = (int32_t*)(tableauBottom.row(y).data);
        uint8_t* errRow = errorBuffer.row(y).data;
        for (int x = 0; x < params.xMax; x++)
        {
            if (params.salientPoints and salientBuffer(y, x) == 0)
            {
                for (int hypIdx = 0; hypIdx < params.hypMax; hypIdx++)
                {
                    smallDisparity(y, x * params.hypMax + hypIdx) = -1;
                }
                continue;
            }
            int minCost = 0;
            int minCostDisp = -1;
            
            
            //compute the cost vector
//            vector<int> costVec;
//            costVec.reserve(params.dispMax);
//            int base = x * params.dispMax;
//            for (int d = 0; d < params.dispMax; d++)
//            {
//                costVec.push_back(dynRow1[base + d] + dynRow2[base + d] 
//                            + dynRow3[base + d] + dynRow4[base + d] - 2*errRow[base + d]);
//            }
            
            
           /* //compute the minima
            vector<int> minIdxVec = findLocalMinima(costVec);
            
            //sort minima
            vector<pair<int, int>> indexedCostVec;
            indexedCostVec.reserve(minIdxVec.size());
            for (auto & idx : minIdxVec)
            {
                indexedCostVec.emplace_back(costVec[idx], idx);
            }
//            sort(indexedCostVec.begin(), indexedCostVec.end());
            sizeAcc += indexedCostVec.size();
            sizeCount++;
            partial_sort(indexedCostVec.begin(), indexedCostVec.begin() + params.hypMax, indexedCostVec.end());
            for (int hypIdx = 0; hypIdx < minIdxVec.size(); hypIdx++)
            {
                smallDisparity(y, x * params.hypMax + hypIdx) = indexedCostVec[hypIdx].second;
                finalErrorMat(y, x * params.hypMax + hypIdx) = indexedCostVec[hypIdx].first;
            }*/
            bool verbose = (y == 0 and x == 24);
            for (int hypIdx = 0; hypIdx < params.hypMax; hypIdx++)
            {
                int32_t & bestDisp = smallDisparity(y, x * params.hypMax + hypIdx);
                int32_t & bestErr = finalErrorMat(y, x * params.hypMax + hypIdx);
                bestErr = INT32_MAX;
                bestDisp = -1;
                int acc1 = -1, acc2 = -1, acc3 = -1;
                for (int d = 0; d < params.dispMax; d++)
                {
                    int base = x * params.dispMax;
                    if (errRow[base + d] > params.maxError) continue;
                    acc1 = acc2;
                    acc2 = acc3;
                    acc3 = dynRow1[base + d] + dynRow2[base + d] 
                            + dynRow3[base + d] + dynRow4[base + d] - 2*errRow[base + d];
                            
                    bool localMin = false;
                    if ( acc2 == -1) continue;
                    else if (acc1 == -1)
                    { 
                        localMin = (acc2 < acc3);
                    }
                    else
                    {
                        localMin = (acc2 < acc3 and acc2 <= acc1);
                    } 
                    if (not localMin) continue;
                    int d2 = d - 1;    
                    if ( bestErr > acc2 and ( acc2 > minCost or 
                            (acc2 == minCost and d2 != minCostDisp) ) )
                    {
                        if (hypIdx > 0 and acc2 > minCost + params.maxHypDiff) continue;
                        bestDisp = d2;
                        bestErr = acc2;
                    }
                }
                if (bestDisp == -1)
                {
                    for (int h2Idx = hypIdx + 1; h2Idx < params.hypMax; h2Idx++)
                    {
                        smallDisparity(y, x * params.hypMax + h2Idx) = -1;
                    } 
                    break;
                }
                minCost = bestErr;
                minCostDisp = bestDisp;
            }
            if (params.verbosity > 4) cout << "    x: " << x << " best error: " << finalErrorMat(y, x) << endl;
        }
        if (params.verbosity > 3) cout << "    y: " << y << endl;
    }
}

//TODO remove this function, depricated
void EnhancedSGM::computeDepth(Mat32f & distance)
{
    if (params.verbosity > 0) cout << "EnhancedSGM::computeDepth(Mat32f &)" << endl;
    distance.create(params.yMax, params.xMax);
    for (int y = 0; y < params.yMax; y++)
    {
        for (int x = 0; x < params.xMax; x++)
        {
            distance(y, x) = computeDepth(x, y);
        }
    }
}

bool EnhancedSGM::computeDepth(double & dist, double & sigma, int x, int y, int h)
{
    if (params.verbosity > 3) 
    {
        cout << "EnhancedSGM::computeDepth(double & dist, double & sigma, int x, int y, int h)" << endl;
    }
    assert(h < params.hypMax);
    
    int idx = getLinearIndex(x, y);
    if (not maskVec[idx])
    { 
        dist = OUT_OF_RANGE;
        sigma = OUT_OF_RANGE;
        return false;
    
    }
    int disparity = smallDisparity(y, x*params.hypMax + h);
    
    if (disparity < 0) 
    {
        dist = OUT_OF_RANGE;
        sigma = OUT_OF_RANGE;
        return true;
    }
    else if (disparity == 0)
    {
        dist = MAX_DEPTH;
//        dist = 1 /MAX_DEPTH; //FIXME
        sigma = MAX_DEPTH / 7; //FIXME
        return true;
    }
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    // to compute point on the second image
    // TODO make virtual rasterizer using the cache
    int u21, v21, u22, v22;
    if (params.useUVCache)
    {
        const int uvCacheStep = params.dispMax + 2 * DISPARITY_MARGIN;
        u21 = uCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity);
        u22 = uCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity + 1);
        v21 = vCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity);
        v22 = vCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity + 1);
    }
    else
    {       
        CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
        raster.steps(disparity);
        u21 = raster.u;
        v21 = raster.v;
        raster.step();
        u22 = raster.u;
        v22 = raster.v;
    }
    
    dist = triangulate(pt1[0], pt1[1], u21, v21);
//    if (dist < 0)
//    {
//        cout << y << " " << x << "    " << dist << " " << disparity << endl;
//        cout << pt1.transpose() << "   " << u21  << " " << v21 << endl;
//    }
    if (dist > MIN_DEPTH)
    {
        double d2 = triangulate(pt1[0], pt1[1], u22, v22);
        
        if (d2 > MIN_DEPTH)
        {
            sigma = abs(d2 - dist) * SIGMA_COEFF;
            return true;
        }
    }
    //if any of triangulations is not computed
    dist = OUT_OF_RANGE;
    sigma = OUT_OF_RANGE;
    return false;
}

double EnhancedSGM::computeDepth(int x, int y, int h)
{
    if (params.verbosity > 3) cout << "EnhancedSGM::computeDepth(int x, int y, int h)" << endl;
    assert(not (params.hypMax == 1 and h > 0));
    
    int idx = getLinearIndex(x, y);
    if (not maskVec[idx]) return 0;
    int disparity = smallDisparity(y, x*params.hypMax + h);
    if (disparity < 0) 
    {
        return OUT_OF_RANGE;
    }
    
    // to compute point on the second image
    int u2, v2;
    if (params.useUVCache)
    {
        const int uvCacheStep = params.dispMax + 2 * DISPARITY_MARGIN;
        u2 = uCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity);
        v2 = vCache(y, x*uvCacheStep + DISPARITY_MARGIN + disparity);
    }
    else
    {    
        CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser2(idx);
        raster.steps(disparity);
        u2 = raster.u;
        v2 = raster.v;
    }
    // point on the first image
    const auto & pt1 = pointVec1[idx];
    
    return triangulate(pt1[0], pt1[1], u2, v2);
}


