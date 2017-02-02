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
Depth-from-motion class for semidense depth estimation
*/


#include "reconstruction/eucm_motion_stereo.h"

#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "utils/filter.h"
#include "geometry/geometry.h"
#include "projection/eucm.h"

#include "reconstruction/stereo_misc.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/epipolar_descriptor.h"


//TODO cam1, cam2 are not needed
void MotionStereo::reprojectDepth(Transf T12, const Mat8u & img2, DepthMap & depth)
{
    
    setTransformation(T12);
    
    const int LENGTH = _params.descLength;
    const int HALF_LENGTH = LENGTH / 2;
    
    
    MHPack flatPack, salientPack;
    vector<vector<uint8_t>> descriptorVec;
    for (int y = 0; y < depth.yMax; y++)
    {
        for (int x = 0; x < depth.xMax; x++)
        {
            Vector2d pt(depth.uConv(x), depth.vConv(y));
            Vector3d X;
            if (not _camera1->reconstructPoint(pt, X)) continue;
            
            CurveRasterizer<int, Polynomial2> descRaster(round(pt), epipoles().getFirstPx(),
                                            _epipolarCurves.getFirst(X));
            if (epipoles().firstIsInverted()) descRaster.setStep(-1);
            vector<uint8_t> descriptor;
            const int step = _epipolarDescriptor.compute(_img1, descRaster, descriptor);
            if (step < 1) continue;
            
            if (not _epipolarDescriptor.goodResp())
            {
                // cannot be used for motion stereo
                
                // the depth is not defined either
                if (depth.at(x, y) < MIN_DEPTH) continue; 
                
                flatPack.imagePointVec.push_back(pt);
            }
            else
            {
//                    if (depth.at(x, y) >= MIN_DEPTH) continue;  //FIXME
                descriptorVec.push_back(descriptor);
                salientPack.imagePointVec.push_back(pt);
            }
        }
    }
    
    //TODO check that sigmaVec is reconstructed
    depth.reconstruct(flatPack, QUERY_POINTS | SIGMA_VALUE);
    depth.reconstruct(salientPack,
            QUERY_POINTS /*| DEFAULT_VALUES*/  | MINMAX | ALL_HYPOTHESES | SIGMA_VALUE | INDEX_MAPPING);
    depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
    
    // for the flat pack project points and replace the hypotheses in depth
    
    T12.inverseTransform(flatPack.cloud, flatPack.cloud);
    cout << flatPack.cloud.size() << "  " << flatPack.sigmaVec.size() << endl;
    for (int idx = 0; idx < flatPack.cloud.size(); idx++)
    {
        depth.pushHypothesis(flatPack.cloud[idx], flatPack.sigmaVec[idx]);
    }
    
    // for the salient pack compute stereo and for corresponding pixel push new hypothesis
    Vector3dVec cloud2;
    T12.inverseTransform(salientPack.cloud, cloud2);
    int dispAcc = 0;
    int dispCount = 0;
    for (int idx = 0; idx < salientPack.imagePointVec.size(); idx++)
    {
        
        // project min-max points
        Vector2d ptMin, ptMax;
        _camera2->projectPoint(cloud2[2*idx], ptMin);
        _camera2->projectPoint(cloud2[2*idx + 1], ptMax);
        
//            cout << cloud2[2*idx].norm() << " " << cloud2[2*idx + 1].norm() << endl;
        
        // if distance is small push depth hyp with the same sigma
        if ((ptMin - ptMax).squaredNorm() < 5)
        {
            depth.pushHypothesis(0.5*(cloud2[2*idx] + cloud2[2*idx + 1]),
                        salientPack.sigmaVec[idx]);
        }
        else    //TODO make a separate function to compute stereo and merge it
                // call it from here and from computeDepth
        {
            // if distance is big enough
            // search along epipolar curve
            //TODO optimize fo Vector2i 
            
            // query 3D point must be projected into 1st frame
            CurveRasterizer<int, Polynomial2> raster(round(ptMax), round(ptMin),
                                            _epipolarCurves.getSecond(salientPack.cloud[2*idx]));
            
            Vector2i diff = round(ptMax - ptMin);
            const int distance = min(int(diff.norm()), _params.dispMax);
            
            raster.steps(-HALF_LENGTH);
            vector<uint8_t> sampleVec;
            vector<int> uVec, vVec;
            const int margin = LENGTH - 1;
            uVec.reserve(distance + margin);
            vVec.reserve(distance + margin);
            sampleVec.reserve(distance + margin);
            for (int d = 0; d < distance + margin; d++, raster.step())
            {
                if (raster.v < 0 or raster.v >= img2.rows 
                    or raster.u < 0 or raster.u >= img2.cols) sampleVec.push_back(0);
                else sampleVec.push_back(img2(raster.v, raster.u));
                uVec.push_back(raster.u);
                vVec.push_back(raster.v);
            }
            
            vector<uint8_t> & descriptor = descriptorVec[salientPack.idxMapVec[idx]];
            
            vector<int> costVec = compareDescriptor(descriptor, sampleVec, _params.flawCost);
            
//                int dBest = 0;
//                int eBest = LENGTH*65535;
//                int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
//                for (int d = 0; d < distance; d++)
//                {
//                    int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                    int bias = 0; //min(_params.biasMax, max(-_params.biasMax, (sum2 - sum1) / NORMALIZER));
//                    int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
//                                        descriptor.begin(), sampleVec.begin() + d, bias, 1);
//                    if (eBest > acc)
//                    {
//                        eBest = acc;
//                        dBest = d;
//                    }
//                }
//                dispAcc += dBest;
//                dispCount++;



                //FIXME
                int dBest = 0;
                int eBest = LENGTH*65535;
                for (int d = 0; d < distance; d++)
                {
                    int acc = costVec[d + HALF_LENGTH];
                    if (eBest > acc)
                    {
                        eBest = acc;
                        dBest = d;
                    }
                }
        
        
        


//                auto poly1 = _epipolarCurves.getSecond(salientPack.cloud[2*idx]);
//                
//                cout << "    curve : " << poly1(ptMin[0], ptMin[1]) << " " 
//                    << poly1(ptMax[0], ptMax[1]) <<  " " 
//<< poly1(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH]) <<endl;
            
//                cout << ptMax.transpose() << " / " << uVec[dBest + HALF_LENGTH] << " " 
//                    << vVec[dBest + HALF_LENGTH] << " / " << ptMin.transpose() 
//<< " / " << dBest << " " << distance << endl;
            
            // put the original depth to the new pixel
       
//                int xd = depth.xConv(uVec[dBest + HALF_LENGTH]);
//                int yd = depth.yConv(vVec[dBest + HALF_LENGTH]);
//                if (not depth.isValid(xd, yd)) continue;
//                depth.at(xd, yd) = 0.5*(salientPack.cloud[2*idx] + salientPack.cloud[2*idx + 1]).norm();
//                depth.sigma(xd, yd) = 1;

            // put the original hypothesis
//                depth.pushHypothesis(0.5*(salientPack.cloud[2*idx] + salientPack.cloud[2*idx + 1]), 1);

            // triangulate and improve sigma
            double d1 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                    uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], CAMERA_2);
            double d2 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                    uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], CAMERA_2);
            double sigma1 = abs(d2 - d1)/2;
            depth.pushImageHypothesis(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], d1, sigma1);
            
//                if ( depth.nearestSigma(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH]) > 0.5)
//                {
//                    cout << distance << endl;
//                    cout << "depth: " << depth.nearest(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH])
//                        << " +-" << depth.nearestSigma(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH]) 
//                        << endl;
//                    cout << "samples:" << endl;
//                    for (auto & x : sampleVec)
//                    {
//                        cout << " " << int(x);
//                    }
//                    cout << endl;
//                    cout << "descriptor:" << endl;
//                    for (auto & x : descriptor)
//                    {
//                        cout << " " << int(x);
//                    }
//                    cout << endl;
//                    cout << "cost:" << endl;
//                    for (auto & x : costVec)
//                    {
//                        cout << " " << int(x);
//                    }
//                    cout << endl<< endl;
//                }
        }
        
    }     
    cout << double(dispAcc) / dispCount << endl;
    //release dynamic objects
}    

//FIXME the same as DepthMap::filter
void filter(double & v1, double & s1, const double v2, const double s2)
{
    double denom = s1 + s2;
    v1 = (v1 * s2 + v2 * s1) / denom;
    s1 = max(s1 * s2 / denom, 0.1);
}

bool MotionStereo::selectPoint(int x, int y)
{
    gu = _params.uConv(x);
    gv = _params.vConv(y);
    flags |= GLB_UV;
    
    //--Check point's saliency
    if (_maskMat(gv, gu) < _params.gradientThresh) return false;
    
    Vector2d pt(gu, gv);
    
    if (not _camera1->reconstructPoint(pt, gX)) return false;
    flags |= GLB_X;
    
    CurveRasterizer<int, Polynomial2> descRaster(round(pt), epipoles().getFirstPx(),
                                _epipolarCurves.getFirst(gX));
    if (epipoles().firstIsInverted()) descRaster.setStep(-1);
    
    gstep = _epipolarDescriptor.compute(_img1, descRaster, gdescriptor);
    if (gstep < 1 or not _epipolarDescriptor.goodResp()) return false;
    flags |= GLB_STEP | GLB_DESCRIPTOR;
    return true;
}


bool MotionStereo::computeUncertainty(double d, double s)
{
    //TODO replace assert?
    uint32_t neededFlag = GLB_X;
    assert( (flags & neededFlag) ^ neededFlag == 0);
    
    Vector3d Xmax;
    
    if (d == OUT_OF_RANGE)  // no prior
    {
        //just rotate
        Xmax = R21() * gX;
        if (not _camera2->projectPoint(Xmax, gptStart)) return false;
        gdispMax = _params.dispMax;
        flags |= GLB_START_POINT | GLB_DISP_MAX;
    }
    else // there is a prior
    {
        gX.normalize();
        Xmax = gX * (d + 2 * s);
        Vector3d Xmin = gX * max(d - 2 * s, MIN_DEPTH);
        Xmax = R21() * (Xmax - t12());
        Xmin = R21() * (Xmin - t12());
        Vector2d ptFin;
        if (not _camera2->projectPoint(Xmax, gptStart)) return false;
        if (not _camera2->projectPoint(Xmin, ptFin)) return false;
        int delta = round((ptFin - gptStart).norm());
        gdispMax = min( _params.dispMax, delta + 1);
        flags |= GLB_START_POINT | GLB_DISP_MAX;
    }
    return true;
}

bool MotionStereo::sampleImage(const Mat8u & img2)
{
    uint32_t neededFlag = GLB_START_POINT | GLB_DISP_MAX | GLB_STEP | GLB_X;
    assert(flags & neededFlag == neededFlag);
    
    Vector2i ptStartRound = round(gptStart);
    if ((ptStartRound - epipoles().getSecondPx()).squaredNorm() < 2500) return false;
    int distance = gdispMax / gstep + MARGIN;
    
    CurveRasterizer<int, Polynomial2> raster(ptStartRound, epipoles().getSecondPx(),
                                _epipolarCurves.getSecond(gX));
    //Important : Epipolar curves are accessed by the reconstructed point in the FIRST frame
                                
    if (epipoles().secondIsInverted()) raster.setStep(-1);
    raster.setStep(gstep);
    raster.steps(-HALF_LENGTH);

    guVec.clear();
    guVec.reserve(distance);
    gvVec.clear();
    gvVec.reserve(distance);
    gsampleVec.clear();
    gsampleVec.reserve(distance);
    for (int d = 0; d < distance; d++, raster.step())
    {
        if (raster.v < 0 or raster.v >= img2.rows 
            or raster.u < 0 or raster.u >= img2.cols) 
        {
            return false;
        }//sampleVec.push_back(0);
        
        gsampleVec.push_back(img2(raster.v, raster.u));
        guVec.push_back(raster.u);
        gvVec.push_back(raster.v);
    }
    flags |= GLB_SAMPLE_VEC | GLB_UV_VEC;
    return true;
}

void MotionStereo::reconstructFirst(double & dist, double & sigma, double & cost)
{
    uint32_t neededFlag = GLB_SAMPLE_VEC | GLB_UV | GLB_UV_VEC | GLB_DESCRIPTOR;
    assert(flags & neededFlag == neededFlag);
    
    vector<int> costVec = compareDescriptor(gdescriptor, gsampleVec, _params.flawCost);
    auto bestCostIter = min_element(costVec.begin(), costVec.end());
    
//    if (gu < 550 and gu > 450 and gv > 280 and gv < 380 )
//    {
//        cout << gu << "   " << gv << endl;
//        cout << "depth: " << dist
//            << " +-" << sigma
//            << endl;
//        cout << "samples:" << endl;
//        for (auto & x : gsampleVec)
//        {
//            cout << " " << int(x);
//        }
//        cout << endl;
//        cout << "coordinates:" << endl;
//        for (auto & x : guVec)
//        {
//            cout << " " << int(x);
//        }
//        cout << endl;
//        for (auto & x : gvVec)
//        {
//            cout << " " << int(x);
//        }
//        cout << endl;
//        cout << "descriptor:" << endl;
//        for (auto & x : gdescriptor)
//        {
//            cout << " " << int(x);
//        }
//        cout << endl;
//        cout << "cost:" << endl;
//        for (auto & x : costVec)
//        {
//            cout << " " << int(x);
//        }
//        cout << endl<< endl;
//    }
    
    if (*bestCostIter < _params.maxError)
    {
        int dBest = bestCostIter - costVec.begin();
        triangulate(gu, gv, guVec[dBest], gvVec[dBest], guVec[dBest + 1], gvVec[dBest + 1],
                dist, sigma, CAMERA_1); 
//        if (sigma > 1 and gu < 600 and gu > 400 and gv > 300 and gv < 550 )
//        {
//            cout << sigma << " " << dist << endl;
//            cout    << gu << " " << gv << " " 
//                    <<  guVec[dBest] << " "  << gvVec[dBest] << " " 
//                     << guVec[dBest + 1] << " "  << gvVec[dBest + 1] << endl;
//        }
        cost = *bestCostIter;
    }
}

void MotionStereo::reconstructSecond(DepthMap & depth)
{
    uint32_t neededFlag = GLB_SAMPLE_VEC | GLB_UV | GLB_UV_VEC | GLB_DESCRIPTOR;
    vector<int> costVec = compareDescriptor(gdescriptor, gsampleVec, _params.flawCost);
    auto bestCostIter = min_element(costVec.begin(), costVec.end());
    
    if (*bestCostIter < _params.maxError)
    {
        int dBest = bestCostIter - costVec.begin();
        int x = depth.xConv(guVec[dBest]);
        int y = depth.yConv(gvVec[dBest]);
        
        double & dist = depth.at(x, y);
        double & sigma = depth.sigma(x, y); 
        double & cost = depth.cost(x, y);
        
        triangulate(gu, gv, guVec[dBest], gvVec[dBest], guVec[dBest + 1], gvVec[dBest + 1],
                dist, sigma, CAMERA_2); 
                //TODO does not work for CAMERA_2 yet
        cost = *bestCostIter;
    }
}

DepthMap MotionStereo::compute(Transf T12, const Mat8u & img2)
{
    //init necessary data structures
    setTransformation(T12);
    DepthMap depthOut(_camera1, _params);
    depthOut.setTo(OUT_OF_RANGE, OUT_OF_RANGE, OUT_OF_RANGE);

    //for each point
    for (int y = 0; y < depthOut.yMax; y++)
    {
        for (int x = 0; x < depthOut.xMax; x++)
        {
            flags = 0;
            if (not selectPoint(x, y)) continue;
            
            if (not computeUncertainty(OUT_OF_RANGE, OUT_OF_RANGE)) continue;
            
            if (not sampleImage(img2)) continue;
            
            reconstructFirst(depthOut.at(x, y), depthOut.sigma(x, y), depthOut.cost(x, y));
        }
    }
    return depthOut;
}

DepthMap MotionStereo::compute(Transf T12, const Mat8u & img2, const DepthMap & depthIn, int arg)
{
    //init necessary data structures
    setTransformation(T12);
    assert(ScaleParameters(depthIn) == ScaleParameters(_params));
    DepthMap depthOut = depthIn;
//    depthOut.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    //for each point
    for (int y = 0; y < depthOut.yMax; y++)
    {
        for (int x = 0; x < depthOut.xMax; x++)
        {
//            depthOut.at(x, y) *= 2;
//            if (arg != 0 and depthIn.at(x, y) != OUT_OF_RANGE) continue;
            flags = 0;
            if (not selectPoint(x, y)) { count1++; continue;}
            
            if (not computeUncertainty(depthIn.at(x, y), depthIn.sigma(x, y))) { count2++; continue;}
            
            if (gdispMax / gstep < 2) 
            {
                count3++;
                depthOut.at(x, y) = depthIn.at(x, y);
                depthOut.sigma(x, y) = depthIn.sigma(x, y);
                continue;  
            }  // the uncertainty is too small
            
            //TODO if the uncertainty is small, fuse the two measurements 
            // replace the old one otherwise
            
            if (not sampleImage(img2)) { count4++; continue;}
            
            reconstructFirst(depthOut.at(x, y), depthOut.sigma(x, y), depthOut.cost(x, y));
        }
    }
    cout << count1 << endl;
    cout << count2 << endl;
    cout << count3 << endl;
    cout << count4 << endl;
    return depthOut;
}
