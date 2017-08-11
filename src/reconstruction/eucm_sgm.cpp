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
#include "projection/eucm.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

CurveRasterizer<int, Polynomial2> EnhancedSgm::getCurveRasteriser(CameraIdx camIdx, 
        int idx, uint32_t * flags) const
{
    Vector2i pti;
    if (camIdx == CAMERA_1) pti = _pointPxVec1[idx];
    else if (camIdx == CAMERA_2) pti = _pinfPxVec[idx]; 
    uint32_t useInverted = epipoles().chooseEpipole(camIdx, pti, _params.epipoleMargin);
    Vector2i goal = epipoles().getPx(camIdx, useInverted);
    CurveRasterizer<int, Polynomial2> raster(pti, goal,
                            _epipolarCurves.get(camIdx, _reconstVec[idx]));
    if (useInverted & EPIPOLE_INVERTED) raster.setStep(-1);
    if (flags != NULL) *flags = useInverted;
    return raster;
}


//TODO reconstruct the depth points, not everything
void EnhancedSgm::computeReconstructed()
{
    _pointVec1.resize(_params.yMax*_params.xMax);
    _pointPxVec1.resize(_params.yMax*_params.xMax);
    for (int y = 0; y < _params.yMax; y++)
    {
        for (int x = 0; x < _params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            _pointVec1[idx] = Vector2d(_params.uConv(x), _params.vConv(y));
            _pointPxVec1[idx] = Vector2i(_params.uConv(x), _params.vConv(y));
        }
    }
    _camera1->reconstructPointCloud(_pointVec1, _reconstVec, _maskVec);
}

void EnhancedSgm::computeRotated()
{
    transf().inverseRotate(_reconstVec, _reconstRotVec);
}

//FIXME _maskVec must be recomputed to discard not projected pInf
void EnhancedSgm::computePinf()
{
    _camera2->projectPointCloud(_reconstRotVec, _pinfVec);
    _pinfPxVec.resize(_pinfVec.size());
    for (int i = 0; i < _pinfVec.size(); i++)
    {
        if (not _maskVec[i]) continue;
        _pinfPxVec[i] = round(_pinfVec[i]);
    }
}

void EnhancedSgm::computeUVCache()
{
    for (int y = 0; y < _params.yMax; y++)
    {
        for (int x = 0; x < _params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            if (not _maskVec[idx]) continue;
            CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(CAMERA_2, idx);
            raster.steps(-DISPARITY_MARGIN);
            const int u_vCacheStep = _params.dispMax + 2 * DISPARITY_MARGIN;
            int32_t * uPtr = (int32_t *)_uCache.row(y).data + x*u_vCacheStep;
            int32_t * vPtr = (int32_t *)_vCache.row(y).data + x*u_vCacheStep;
            for (int i = 0; i  < u_vCacheStep; i++, raster.step(), uPtr++, vPtr++)
            {
                if (raster.v < 0 or raster.v >= _params.vMax 
                    or raster.u < 0 or raster.u >= _params.uMax)
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

void EnhancedSgm::createBuffer()
{
    if (_params.verbosity > 1) cout << "EnhancedSgm::createBuffer" << endl;
    assert(_params.hypMax > 0);
    int bufferWidth = _params.xMax*_params.dispMax;
    _stepBuffer.create(_params.yMax, _params.xMax);
    _errorBuffer.create(_params.yMax, bufferWidth);
    _tableauLeft.create(_params.yMax, bufferWidth);
    _tableauRight.create(_params.yMax, bufferWidth);
    _tableauTop.create(_params.yMax, bufferWidth);
    _tableauBottom.create(_params.yMax, bufferWidth);
    _smallDisparity.create(_params.yMax, _params.xMax * _params.hypMax);
    _finalErrorMat.create(_params.yMax, _params.xMax * _params.hypMax);
    _skipBuffer.create(_params.yMax, _params.xMax);
    if (_params.imageBasedCost) _costBuffer.create(_params.yMax, _params.xMax);
    if (_params.salientPoints) _salientBuffer.create(_params.yMax, _params.xMax);
    _uCache.create(_params.yMax, _params.xMax * (_params.dispMax + 2*DISPARITY_MARGIN));
    _vCache.create(_params.yMax, _params.xMax * (_params.dispMax + 2*DISPARITY_MARGIN));
    if (_params.verbosity > 2) 
    {
        cout << "    small disparity size: " << _smallDisparity.size() << endl;
    }
}

void EnhancedSgm::computeStereo(const Mat8u & img1, const Mat8u & img2, DepthMap & depth)
{
    _skipBuffer.setTo(0);
    computeCurveCost(img1, img2);
    
    computeDynamicProgramming();
    
    if (_params.hypMax == 1) reconstructDisparity();
    else reconstructDisparityMH();
    
    reconstructDepth(depth);    

}

void EnhancedSgm::reconstructDepth(DepthMap & depth) const
{
    if (_params.verbosity > 2) 
    {
        cout << "EnhancedSgm::reconstructDepth(DepthMap & depth)" << endl;
    }
    depth = DepthMap(_camera1, _params, _params.hypMax);
    for (int h = 0; h < _params.hypMax; h++)
    {
        for (int y = 0; y < _params.yMax; y++)
        {
            for (int x = 0; x < _params.xMax; x++)
            {
                if (_params.salientPoints and not _salientBuffer(y, x) or _skipBuffer(y, x)) 
                {
                    depth.at(x, y, h) = OUT_OF_RANGE;
                    depth.sigma(x, y, h) = OUT_OF_RANGE;
                    depth.cost(x, y, h) = OUT_OF_RANGE;
                    continue;
                }
                depth.cost(x, y, h) = _errorBuffer(y, x*_params.dispMax + h);
                
                int idx = getLinearIndex(x, y);
                if (not _maskVec[idx])
                { 
                    depth.at(x, y, h) = OUT_OF_RANGE;
                    depth.sigma(x, y, h) = OUT_OF_RANGE;
                    depth.cost(x, y, h) = OUT_OF_RANGE;
                    continue;
                }
                int disparity = _smallDisparity(y, x*_params.hypMax + h);
                
                // point on the first image
                const auto & pt1 = _pointVec1[idx];
                
                // to compute point on the second image
                // TODO make virtual rasterizer using the cache
                int u21, v21, u22, v22;
                int step = _stepBuffer(y, x);
                if (_params.useUVCache)
                {
                    const int u_vCacheStep = _params.dispMax + 2 * DISPARITY_MARGIN;
                    u21 = _uCache(y, x*u_vCacheStep + DISPARITY_MARGIN + disparity);
                    u22 = _uCache(y, x*u_vCacheStep + DISPARITY_MARGIN + disparity + step);
                    v21 = _vCache(y, x*u_vCacheStep + DISPARITY_MARGIN + disparity);
                    v22 = _vCache(y, x*u_vCacheStep + DISPARITY_MARGIN + disparity + step);
                }
                else
                {       
                    CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(CAMERA_2, idx);
                    raster.steps(disparity);
                    u21 = raster.u;
                    v21 = raster.v;
                    raster.steps(step);
                    u22 = raster.u;
                    v22 = raster.v;
                }
                
                triangulate(pt1[0], pt1[1], u21, v21, u22, v22,
                             depth.at(x, y, h), depth.sigma(x, y, h));
                        
            }
        }
    }
}

void EnhancedSgm::skipPixel(int x, int y)
{
    uint8_t * outPtr = _errorBuffer.row(y).data + x*_params.dispMax;            
    *outPtr = 0;
    _skipBuffer(y, x) = 1;
    fill(outPtr + 1, outPtr + _params.dispMax, 255);
}

void EnhancedSgm::computeCurveCost(const Mat8u & img1, const Mat8u & img2)
{
    if (_params.verbosity > 0) cout << "EnhancedSgm::computeCurveCost" << endl;
    
    // compute the weights for matching cost
    
    if (_params.salientPoints) _salientBuffer.setTo(0);
    
    for (int y = 0; y < _params.yMax; y++)
    {
        for (int x = 0; x < _params.xMax; x++)
        {
            int idx = getLinearIndex(x, y);
            if (_params.verbosity > 5) 
            {
                cout << "    x: " << x << " y: " << y << "  idx: " << idx; 
                cout << "  mask: " << _maskVec[idx] <<  endl;
            }
            if (not _maskVec[idx])
            {
                skipPixel(x, y);
                continue;
            }
            // compute the local image descriptor,
            // a piece of the epipolar curve on the first image
            vector<uint8_t> descriptor;
            uint32_t flags;
            CurveRasterizer<int, Polynomial2> descRaster = getCurveRasteriser(CAMERA_1, idx, &flags);
            if (flags & EPIPOLE_TOO_CLOSE) 
            {
                skipPixel(x, y);
                continue;
            }
            const int step = _epipolarDescriptor.compute(img1, descRaster, descriptor);
            _stepBuffer(y, x) = step;
            if (step < 1) 
            {
                skipPixel(x, y);
                continue;
            }
            if (_params.imageBasedCost) 
            {
                switch (step)
                {
                case 1:
                    _costBuffer(y, x) = _params.lambdaJump;
                    break;
                case 2:
                    _costBuffer(y, x) = _params.lambdaJump * 3;
                    break;
                default:
                    _costBuffer(y, x) = _params.lambdaJump * 6;
                    break;
                }
            }
            
            //TODO revise the criterion (step == 1)
            if (_params.salientPoints and step <= 2 and _epipolarDescriptor.goodResp())
            {
                _salientBuffer(y, x) = 1;
            }
            const int nSteps = ( _params.dispMax  + step - 1 ) / step; 
               
            //sample the curve 
            vector<uint8_t> sampleVec(nSteps + MARGIN, 0);
            bool crossedImageBoundary = false;
            if (_params.useUVCache)
            {
                const int u_vCacheStep = _params.dispMax + 2 * DISPARITY_MARGIN;
                int32_t * uPtr = (int32_t *)_uCache.row(y).data + x*u_vCacheStep;
                int32_t * vPtr = (int32_t *)_vCache.row(y).data + x*u_vCacheStep;
                uPtr += DISPARITY_MARGIN - HALF_LENGTH * step;
                vPtr += DISPARITY_MARGIN - HALF_LENGTH * step;
                for (int i = 0; i  < nSteps + MARGIN; i++, uPtr += step, vPtr += step)
                {
                    if (*uPtr < 0 or *vPtr < 0) 
                    {
                        crossedImageBoundary = true;
                        break;
                    }
                    else sampleVec[i] = img2(*vPtr, *uPtr);
                }
            }
            else
            {
                CurveRasterizer<int, Polynomial2> raster = getCurveRasteriser(CAMERA_2, idx);
                raster.setStep(step); 
                raster.steps(-HALF_LENGTH);           
                
                if (_params.verbosity > 6)
                {
                    cout << "CURVE RASTERIZER" << endl;
                    cout << "delta : " << raster.delta << endl;
                    cout << "fu fv : " << raster.fu << " " << raster.fv << endl;
                    cout << " u  v : " << raster.u << " " << raster.v << endl;
                    
                    const auto & surf = raster.surf;
                    cout << " SURF : " << endl;
                    cout << surf.kuu << " " << surf.kuv << " " << surf.kvv << " " << surf.ku
                         << " " << surf.kv << " " << surf.k1 << endl;
                }
                
                for (int i = 0; i  < nSteps + MARGIN; i++, raster.step())
                {
                    if (raster.v < 0 or raster.v >= img2.rows 
                        or raster.u < 0 or raster.u >= img2.cols)
                        {
                            crossedImageBoundary = true;
                            break;
                        }
                    if (_params.verbosity > 5)
                    {
                        cout << raster.u << "  " << raster.v << endl;
                    }
                    sampleVec[i] = img2(raster.v, raster.u);
                }
            }
            if (crossedImageBoundary)
            {
                skipPixel(x, y);
                continue;
            }
            vector<int> costVec = compareDescriptor(descriptor, sampleVec, _params.flawCost);
            
            if (_params.verbosity > 4)
            {
                cout << "Point : " << x << " " << y << endl;
                cout << "Step : " << step << endl;
                cout << "samples :" << endl;
                for (auto & x : sampleVec)
                {
                    cout << setw(6) << int(x);
                }
                cout << endl;
                cout << "cost :" << endl;
                for (auto & x : costVec)
                {
                    cout << setw(6) << int(x);
                }
                cout << endl;
                cout << "descriptor :" << endl;
                for (auto & x : descriptor)
                {
                    cout << setw(6) << int(x);
                }
                cout << endl;
            }
//            //compute the bias;
//            int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
            
            // fill up the cost buffer
            uint8_t * outPtr = _errorBuffer.row(y).data + x*_params.dispMax;
            auto costIter = costVec.begin() + HALF_LENGTH;
            for (int d = 0; d < nSteps; d++, outPtr += step)
            {
//                int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                int bias = min(_params.maxBias, max(-_params.maxBias, (sum2 - sum1) / LENGTH));
//                int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
//                                descriptor.begin(), sampleVec.begin() + d, bias);
//                *outPtr = acc / NORMALIZER;

                *outPtr = min(*costIter, 255);
                ++costIter;
            }
            if (step > 1) fillGaps(_errorBuffer.row(y).data + x*_params.dispMax, step);
        }
    }
//    cout << "Epipoles : " << endl;
//    cout << epipoles().useInvertedEpipoleSecond(Vector2i(250, 250)) << endl;
//    cout << epipoles().getSecond(true).transpose() << "   
//  "  << epipoles().getSecond(false).transpose() << endl;
//    cout << epipoles().getSecondPx(true).transpose() << "    
// "  << epipoles().getSecondPx(false).transpose() << endl;
//    
//    
//    cout << epipoles().epipole1projected<< epipoles().epipole2projected
//        << epipoles().antiEpipole1projected<< epipoles().antiEpipole2projected << endl;
}

void EnhancedSgm::fillGaps(uint8_t * const data, const int step)
{
    assert(step > 0);
    //linear interpolation for all intermediate points
    int base;
    switch (step)
    {
    case 2:
        for (base = 2; base < _params.dispMax; base += 2)
        {
            data[base - 1] = (data[base - 2] + data[base]) / 2;
        }
        break;
    case 3:
        for (base = 3; base < _params.dispMax; base += 3)
        {
            const uint8_t & val1 = data[base - 3];
            const uint8_t & val2 = data[base];
            data[base - 2] = (val1 << 1  + val2) / 3;
            data[base - 1] = (val1 + val2 << 1) / 3;
        }
        break;
    default:
        for (base = step; base < _params.dispMax; base += step)
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
    for (int i = base + 1 ; i < _params.dispMax; i++)
    {
        data[i] = val;
    }
}

void EnhancedSgm::computeDynamicStep(const int32_t* inCost, const uint8_t * error, int32_t * outCost)
{
    int bestCost = inCost[0];
    for (int i = 1; i < _params.dispMax; i++)
    {
        bestCost = min(bestCost, inCost[i]);
    }
    int & val0 = outCost[0];
    val0 = inCost[0];
    val0 = min(val0, inCost[1] + _params.lambdaStep);
    val0 = min(val0, bestCost + _jumpCost);
    val0 += error[0];
    for (int i = 1; i < _params.dispMax-1; i++)
    {
        int & val = outCost[i];
        val = inCost[i];
        val = min(val, inCost[i + 1] + _params.lambdaStep);
        val = min(val, inCost[i - 1] + _params.lambdaStep);
        val = min(val, bestCost + _jumpCost);
        val += error[i];
    }
    int & vald = outCost[_params.dispMax - 1];
    vald = inCost[_params.dispMax - 1];
    vald = min(vald, inCost[_params.dispMax - 2] + _params.lambdaStep);
    vald = min(vald, bestCost + _jumpCost);
    vald += error[_params.dispMax - 1];
}

void EnhancedSgm::computeDynamicProgramming()
{
    if (_params.verbosity > 0) cout << "EnhancedSgm::computeDynamicProgramming" << endl;
    if (_params.verbosity > 1) cout << "    left" << endl;
    
    if (not _params.imageBasedCost) _jumpCost = _params.lambdaJump;
    
    // left _tableau init
    for (int y = 0; y < _params.yMax; y++)
    {
        int32_t * _tableauRow = (int32_t *)(_tableauLeft.row(y).data);
        uint8_t * errorRow = _errorBuffer.row(y).data;
        // init the first row
        copy(errorRow, errorRow + _params.dispMax, _tableauRow);
        // fill up the _tableau
        for (int x = 1; x < _params.xMax; x++)
        {
            if (_params.imageBasedCost) _jumpCost = _costBuffer(y, x);
            computeDynamicStep(_tableauRow + (x - 1)*_params.dispMax,
                    errorRow + x*_params.dispMax, _tableauRow + x*_params.dispMax);
        }
    }
    if (_params.verbosity > 1) cout << "    right" << endl;  
    // right _tableau init
    for (int y = 0; y < _params.yMax; y++)
    {
        int32_t * _tableauRow = (int32_t *)(_tableauRight.row(y).data);
        uint8_t * errorRow = _errorBuffer.row(y).data;
        int base = (_params.xMax - 1) * _params.dispMax;
        copy(errorRow + base, errorRow + base + _params.dispMax, _tableauRow + base);
        
        for (int x = _params.xMax - 2; x >= 0; x--)
        {
            if (_params.imageBasedCost) _jumpCost = _costBuffer(y, x);
            computeDynamicStep(_tableauRow + (x + 1)*_params.dispMax, 
                    errorRow + x*_params.dispMax, _tableauRow + x*_params.dispMax);
        }
    }
    if (_params.verbosity > 1) cout << "    top" << endl;
    // top-down _tableau init
    for (int x = 0; x < _params.xMax; x++)
    {
        auto _tableauCol = _tableauTop(Rect(x*_params.dispMax, 0, _params.dispMax, _params.yMax));
        auto errorCol = _errorBuffer(Rect(x*_params.dispMax, 0, _params.dispMax, _params.yMax));
        copy(errorCol.data, errorCol.data + _params.dispMax, (int*)(_tableauCol.data));
        for (int y = 1; y < _params.yMax; y++)
        {
            if (_params.imageBasedCost) _jumpCost = _costBuffer(y, x);
            computeDynamicStep((int32_t*)(_tableauCol.row(y-1).data), 
                    errorCol.row(y).data,
                    (int32_t*)(_tableauCol.row(y).data));
        }
    }
    if (_params.verbosity > 1) cout << "    bottom" << endl;
    // bottom-up _tableau init
    for (int x = 0; x < _params.xMax; x++)
    {
        auto _tableauCol = _tableauBottom(Rect(x*_params.dispMax, 0, _params.dispMax, _params.yMax));
        auto errorCol = _errorBuffer(Rect(x*_params.dispMax, 0, _params.dispMax, _params.yMax));
        int vLast = _params.yMax - 1;
        copy(errorCol.row(vLast).data, 
                errorCol.row(vLast).data + _params.dispMax, 
                (int32_t*)(_tableauCol.row(vLast).data));
        for (int y = _params.yMax - 2; y >= 0; y--)
        {
            if (_params.imageBasedCost) _jumpCost = _costBuffer(y, x);
            computeDynamicStep((int32_t*)(_tableauCol.row(y+1).data), 
                    errorCol.row(y).data,
                    (int32_t*)(_tableauCol.row(y).data));
        }
    }
    
}

void EnhancedSgm::reconstructDisparity()
{
    if (_params.verbosity > 0) cout << "EnhancedSgm::reconstructDisparity" << endl;
//    int sizeAcc = 0;
//    int sizeCount = 0;
    for (int y = 0; y < _params.yMax; y++)
    {
        int32_t* dynRow1 = (int32_t*)(_tableauLeft.row(y).data);
        int32_t* dynRow2 = (int32_t*)(_tableauRight.row(y).data);
        int32_t* dynRow3 = (int32_t*)(_tableauTop.row(y).data);
        int32_t* dynRow4 = (int32_t*)(_tableauBottom.row(y).data);
        uint8_t* errRow = _errorBuffer.row(y).data;
        uint8_t* skipRow = _skipBuffer.row(y).data;
        for (int x = 0; x < _params.xMax; x++)
        {
            if ((_params.salientPoints and _salientBuffer(y, x) == 0) or *(skipRow + x))
            {
                _smallDisparity(y, x) = -1;
                continue;
            }
            int32_t & bestDisp = _smallDisparity(y, x);
            int32_t & bestCost = _finalErrorMat(y, x);
            bestCost = INT32_MAX;
            bestDisp = -1;
            const int base = x * _params.dispMax;
            if (_params.verbosity > 4) cout << "Err : ";
            for (int d = 0; d < _params.dispMax; d++)
            {
                const int & err = errRow[base + d];
                if (_params.verbosity > 4) cout << setw(8) << err;
                if (err > _params.maxError) continue;
                int cost = dynRow1[base + d] + dynRow2[base + d] 
                        + dynRow3[base + d] + dynRow4[base + d] - 2 * err;
                
                if ( bestCost > cost)
                {
                    bestDisp = d;
                    bestCost = cost;
                }
            }
             if (_params.verbosity > 4) cout << endl;
            if (_params.verbosity > 4) cout << "    x: " << x << " best error: " 
                    << _finalErrorMat(y, x) << "   d : " << bestDisp <<  endl;
        }
        if (_params.verbosity > 3) cout << "    y: " << y << endl;
    }
}

void EnhancedSgm::reconstructDisparityMH()
{
    if (_params.verbosity > 0) cout << "EnhancedSgm::reconstructDisparityMH" << endl;
    const int hypShift = _params.xMax*_params.yMax;
//    int sizeAcc = 0;
//    int sizeCount = 0;
    for (int y = 0; y < _params.yMax; y++)
    {
        int32_t* dynRow1 = (int32_t*)(_tableauLeft.row(y).data);
        int32_t* dynRow2 = (int32_t*)(_tableauRight.row(y).data);
        int32_t* dynRow3 = (int32_t*)(_tableauTop.row(y).data);
        int32_t* dynRow4 = (int32_t*)(_tableauBottom.row(y).data);
        uint8_t* errRow = _errorBuffer.row(y).data;
        for (int x = 0; x < _params.xMax; x++)
        {
            if (_params.salientPoints and _salientBuffer(y, x) == 0)
            {
                for (int hypIdx = 0; hypIdx < _params.hypMax; hypIdx++)
                {
                    _smallDisparity(y, x * _params.hypMax + hypIdx) = -1;
                }
                continue;
            }
            int minCost = 0;
            int minCostDisp = -1;
            
            
            //compute the cost vector
//            vector<int> costVec;
//            costVec.reserve(_params.dispMax);
//            int base = x * _params.dispMax;
//            for (int d = 0; d < _params.dispMax; d++)
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
            partial_sort(indexedCostVec.begin(), indexedCostVec.begin() + 
                    _params.hypMax, indexedCostVec.end());
            for (int hypIdx = 0; hypIdx < minIdxVec.size(); hypIdx++)
            {
                _smallDisparity(y, x * _params.hypMax + hypIdx) = indexedCostVec[hypIdx].second;
                _finalErrorMat(y, x * _params.hypMax + hypIdx) = indexedCostVec[hypIdx].first;
            }*/
            for (int hypIdx = 0; hypIdx < _params.hypMax; hypIdx++)
            {
                int32_t & bestDisp = _smallDisparity(y, x * _params.hypMax + hypIdx);
                int32_t & bestErr = _finalErrorMat(y, x * _params.hypMax + hypIdx);
                bestErr = INT32_MAX;
                bestDisp = -1;
                int acc1 = -1, acc2 = -1, acc3 = -1;
                for (int d = 0; d < _params.dispMax; d++)
                {
                    int base = x * _params.dispMax;
                    if (errRow[base + d] > _params.maxError) continue;
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
                        if (hypIdx > 0 and acc2 > minCost + _params.maxHypDiff) continue;
                        bestDisp = d2;
                        bestErr = acc2;
                    }
                }
                if (bestDisp == -1)
                {
                    for (int h2Idx = hypIdx + 1; h2Idx < _params.hypMax; h2Idx++)
                    {
                        _smallDisparity(y, x * _params.hypMax + h2Idx) = -1;
                    } 
                    break;
                }
                minCost = bestErr;
                minCostDisp = bestDisp;
            }
            if (_params.verbosity > 4) cout << "    x: " << x << " best error: " 
                << _finalErrorMat(y, x) << endl;
        }
        if (_params.verbosity > 3) cout << "    y: " << y << endl;
    }
}


