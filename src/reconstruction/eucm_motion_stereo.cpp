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
#include "reconstruction/epipoles.h"
#include "reconstruction/eucm_stereo.h"


//TODO cam1, cam2 are not needed
void MotionStereo::reprojectDepth(Transformation<double> T12, const Mat8u & img2, DepthMap & depth)
{
    
    epipolarPtr = new EnhancedEpipolar(T12, camera1, camera2, 2000, params.verbosity);
    StereoEpipoles epipoles(camera1, camera2, T12);
    
    const int LENGTH = params.descLength;
    const int HALF_LENGTH = LENGTH / 2;
    
    vector<int> kernelVec, waveVec;
    const int NORMALIZER = initKernel(kernelVec, LENGTH);
    const int WAVE_NORM = initWave(waveVec, LENGTH);
    EpipolarDescriptor epipolarDescriptor(LENGTH, WAVE_NORM, waveVec.data(), {1});
    
    MHPack flatPack, salientPack;
    vector<vector<uint8_t>> descriptorVec;
    for (int y = 0; y < depth.yMax; y++)
    {
        for (int x = 0; x < depth.xMax; x++)
        {
            Vector2d pt(depth.uConv(x), depth.vConv(y));
            Vector3d X;
            if (not camera1->reconstructPoint(pt, X)) continue;
            
            CurveRasterizer<int, Polynomial2> descRaster(round(pt), epipoles.getFirstPx(),
                                            epipolarPtr->getFirst(X));
            if (epipoles.firstIsInverted()) descRaster.setStep(-1);
            vector<uint8_t> descriptor;
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            if (step < 1) continue;
            
            if (not epipolarDescriptor.goodResp())
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
    Transform12 = T12;
    int dispAcc = 0;
    int dispCount = 0;
    for (int idx = 0; idx < salientPack.imagePointVec.size(); idx++)
    {
        
        // project min-max points
        Vector2d ptMin, ptMax;
        camera2->projectPoint(cloud2[2*idx], ptMin);
        camera2->projectPoint(cloud2[2*idx + 1], ptMax);
        
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
                                            epipolarPtr->getSecond(salientPack.cloud[2*idx]));
            
            Vector2i diff = round(ptMax - ptMin);
            const int distance = min(int(diff.norm()), params.dispMax);
            
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
            
            vector<int> costVec = compareDescriptor(descriptor, sampleVec);
            
//                int dBest = 0;
//                int eBest = LENGTH*65535;
//                int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
//                for (int d = 0; d < distance; d++)
//                {
//                    int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                    int bias = 0; //min(params.biasMax, max(-params.biasMax, (sum2 - sum1) / NORMALIZER));
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
        
        
        


//                auto poly1 = epipolarPtr->getSecond(salientPack.cloud[2*idx]);
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
    delete epipolarPtr;
    epipolarPtr = NULL;
}    

void MotionStereo::computeDepth(Transformation<double> T12, const Mat8u & img2, DepthMap & depth)
{
    if (params.verbosity > 0) cout << "MotionStereo::computeDepth" << endl;
    epipolarPtr = new EnhancedEpipolar(T12, camera1, camera2, 2000, params.verbosity);
    StereoEpipoles epipoles(camera1, camera2, T12);
    Transform12 = T12;
    
    if (params.verbosity > 1) cout << "    descriptor kernel selection" << endl;
    const int LENGTH = params.descLength;
    const int HALF_LENGTH = LENGTH / 2;
    
    // compute the weights for matching cost
    vector<int> kernelVec, waveVec;
    const int NORMALIZER = initKernel(kernelVec, LENGTH);
    const int WAVE_NORM = initWave(waveVec, LENGTH);
    EpipolarDescriptor epipolarDescriptor(LENGTH, LENGTH * 3, waveVec.data(), {1, 2, 3});
    
    MHPack salientPack;
    //TODO to optimize make a continuous vector<uint8_t>
    vector<vector<uint8_t>> descriptorVec;
    Vector2iVec depthPointVec;
    vector<int> stepVec;
    if (params.verbosity > 1) cout << "    beginning loop" << endl;
    
    for (int y = 0; y < depth.yMax; y++)
    {
        for (int x = 0; x < depth.xMax; x++)
        {
            Vector2d pt(depth.uConv(x), depth.vConv(y));
            Vector3d X;
            if (not camera1->reconstructPoint(pt, X)) continue;
            
            CurveRasterizer<int, Polynomial2> descRaster(round(pt), epipoles.getFirstPx(),
                                            epipolarPtr->getFirst(X));
            if (epipoles.firstIsInverted()) descRaster.setStep(-1);
            vector<uint8_t> descriptor;
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            if (step < 1) continue;
            
            if (not epipolarDescriptor.goodResp()) continue;
            descriptorVec.push_back(descriptor);
            salientPack.imagePointVec.push_back(pt);
            stepVec.push_back(step);
            depthPointVec.emplace_back(x, y);
        }
    }
    
    if (params.verbosity > 1) cout << "    reconstructing salient points" << endl;
    depth.reconstruct(salientPack,
         QUERY_POINTS | MINMAX | DEFAULT_VALUES | ALL_HYPOTHESES | SIGMA_VALUE | INDEX_MAPPING); 
    
    //TODO make an option of emptying nonsalient points in the depth map
    depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
    
    Vector3dVec cloud2;
    T12.inverseTransform(salientPack.cloud, cloud2);
    
//        int count = 0;
//        double diffAbs = 0;
//        double diffRel = 0;
    
    if (params.verbosity > 1) cout << "    searching for correspondences" << endl;
    for (int idx = 0; idx < salientPack.imagePointVec.size(); idx++)
    {
        
        // project min-max points
        Vector2d ptMin, ptMax;
        //FIXME
        
        if (not camera2->projectPoint(cloud2[2*idx + 1], ptMax)) continue;
        
        int diffLength = 40; //FIXME
        if (camera2->projectPoint(cloud2[2*idx], ptMin))
        {
//            cout << (ptMin - ptMax).transpose() << " ##### ";
        
            Vector2i diff = round(ptMax - ptMin);
            diffLength = max(abs(diff[0]), abs(diff[1])) + 1;
        }
        const int step = stepVec[salientPack.idxMapVec[idx]];
        const auto & ptEpipole2 = epipoles.getSecondPx();
        const auto ptMaxRound = round(ptMax);
        if (diffLength > 2*step and (ptMaxRound - ptEpipole2).squaredNorm() > 2500) 
        {
            // if distance is big enough
            // search along epipolar curve
            
            // query 3D point must be projected into 1st frame
            CurveRasterizer<int, Polynomial2> raster(ptMaxRound, ptEpipole2,
                                            epipolarPtr->getSecond(salientPack.cloud[2*idx]));
            if (epipoles.secondIsInverted()) raster.setStep(-1);
            
            const int distance = min(diffLength, params.dispMax) / step;
            raster.setStep(step);
            raster.steps(-HALF_LENGTH);
            vector<uint8_t> sampleVec;
            vector<int> uVec, vVec;
            const int margin = LENGTH - 1;
            uVec.reserve(distance + margin);
            vVec.reserve(distance + margin);
            sampleVec.reserve(distance + margin);
            bool imageBorder = false;
            for (int d = 0; d < distance + margin; d++, raster.step())
            {
                if (raster.v < 0 or raster.v >= img2.rows 
                    or raster.u < 0 or raster.u >= img2.cols) 
                    {
                        imageBorder = true;
                        break;
                    }//sampleVec.push_back(0);
                else sampleVec.push_back(img2(raster.v, raster.u));
                uVec.push_back(raster.u);
                vVec.push_back(raster.v);
            }
            if (imageBorder) continue;
            vector<uint8_t> & descriptor = descriptorVec[salientPack.idxMapVec[idx]];
            
            //FIXME tmp
//                double tmpAbs = 0, tmpRel = 0;
//                for (int i = 0; i < descriptor.size(); i++)
//                {
//                    tmpAbs += int(descriptor[i]) - int(sampleVec[i]);
//                    tmpRel += descriptor[i] / double(sampleVec[i] + 1);
//                    if (sampleVec[i] == 0)
//                    {
//                        cout << "NULL : " << uVec[i] << "  " << vVec[i] << endl;
//                    }
//                    if (uVec[i] < 0 or vVec[i] < 0)
//                    {
//                        cout << "out of image : " << round(ptMax).transpose() << endl;
//                    }
//                }
//                if (tmpRel > LENGTH * 5) cout << setw(10)  << int(descriptor[HALF_LENGTH]) << setw(10) 
//                        << int(sampleVec[HALF_LENGTH]) << setw(10)
//                        << tmpAbs/ LENGTH << setw(10) << tmpRel/ LENGTH << endl;
//                diffAbs += tmpAbs / LENGTH;
//                diffRel += tmpRel / LENGTH;
//                count++;
//                continue;
            
            vector<int> costVec = compareDescriptor(descriptor, sampleVec);
            
            if (params.verbosity > 2)
            {
                cout << "Point : " << salientPack.imagePointVec[idx].transpose() << endl;
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
            //TODO make it possible to detect multiple hypotheses if there is no prior
            //TODO make this a parameter
            int dBest = -1;
            int eBest = LENGTH*5;
            for (int d = 0; d < distance; d++)
            {
                const int & acc = costVec[d + HALF_LENGTH];
                if (eBest > acc)
                {
                    eBest = acc;
                    dBest = d;
                }
            }
            if (dBest == -1)
            {
                Vector2i depthPt = depthPointVec[salientPack.idxMapVec[idx]];
                //discard the hypothesis
                depth.at(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = OUT_OF_RANGE;
                depth.sigma(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = OUT_OF_RANGE;
                continue;
            }
            // triangulate and improve sigma
            double d1 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                    uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], CAMERA_1);
            double d2 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                    uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], CAMERA_1);
            double sigma1 = abs(d2 - d1) * SIGMA_COEFF; //variance of a uniform distribution
            
            //TODO make push by idx
            Vector2i depthPt = depthPointVec[salientPack.idxMapVec[idx]];
            
//                DepthMap::filter(depth.at(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]), 
//                        depth.sigma(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]), d1, sigma1);
            
            //TODO make the choise of strategy parametric 
            // and divide this function into blocks
            depth.at(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = d1;
            depth.sigma(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = sigma1;
        }
        
    }
//        cout << "diffAbs : " << diffAbs  / count << endl;
//                cout << "diffRel : " << diffRel / count << endl;
    depth.regularize();
}

//FIXME the same as DepthMap::filter
void filter(double & v1, double & s1, const double v2, const double s2)
{
    double denom = s1 + s2;
    v1 = (v1 * s2 + v2 * s1) / denom;
    s1 = max(s1 * s2 / denom, 0.1);
}

void MotionStereo::validateDepth(Transformation<double> T12, const Mat8u & img2, DepthMap & depth)
{
    //constants
    const int LENGTH = params.descLength;
    const int HALF_LENGTH = LENGTH / 2;
    const int MARGIN = LENGTH - 1;
    
    //init necessary data structures
    EpipolarDescriptor epipolarDescriptor(params.descLength, 5, {1, 2, 3, 4});
    StereoEpipoles epipoles(camera1, camera2, T12);
    epipolarPtr = new EnhancedEpipolar(T12, camera1, camera2, 2000, params.verbosity);
    Transform12 = T12;
    
    //split points int toValidate and toCompute
    MHPack toValidatePack;
    
    if (depth.empty()) //init depth
    {
        depth = DepthMap(camera1, params);
        depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE, OUT_OF_RANGE);
    }
    
    //cache toCompute descriptors
    vector<vector<uint8_t>> descriptorVec;
    vector<int> stepVec;
    Vector2iVec depthValidPtVec, depthComputePtVec;
    Vector2dVec toComputePointVec;
    for (int y = 0; y < depth.yMax; y++)
    {
        for (int x = 0; x < depth.xMax; x++)
        {
            Vector2d pt(depth.uConv(x), depth.vConv(y));
            
            if (depth.at(x, y) != OUT_OF_RANGE) 
            {
                toValidatePack.imagePointVec.push_back(pt);
                depthValidPtVec.emplace_back(x, y);
                continue;
            }
            
            Vector3d X;
            if (not camera1->reconstructPoint(pt, X)) continue;
            
            CurveRasterizer<int, Polynomial2> descRaster(round(pt), epipoles.getFirstPx(),
                                            epipolarPtr->getFirst(X));
            if (epipoles.firstIsInverted()) descRaster.setStep(-1);
            vector<uint8_t> descriptor;
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            if (step < 1 or not epipolarDescriptor.goodResp()) continue;
            if (params.verbosity > 1) 
            {
                cout << x << " " << y << " " << epipolarDescriptor.getResp() << endl;
                for (auto & x : descriptor)
                {
                    cout << int(x) << " ";
                }
                cout << endl;
            }
            descriptorVec.push_back(descriptor);
            toComputePointVec.push_back(pt);
            stepVec.push_back(step);
            depthComputePtVec.emplace_back(x, y);
        }
    }
    
    ////for toValidate
    
    //reconstruct with uncertainty
    depth.reconstruct(toValidatePack, QUERY_POINTS | MINMAX | INDEX_MAPPING);
    Vector3dVec toValidateCloud2;
    T12.inverseTransform(toValidatePack.cloud, toValidateCloud2);
        
    for (int idx = 0; idx < toValidatePack.imagePointVec.size(); idx++)
    {
        // project min-max points
        Vector2d ptMin, ptMax;
        
        if (not camera2->projectPoint(toValidateCloud2[2*idx + 1], ptMax)) continue;
        if (not camera2->projectPoint(toValidateCloud2[2*idx], ptMin)) continue;
        
        const auto ptEpipole2 = epipoles.getSecondPx();
        const auto ptMaxRound = round(ptMax);
        
        //the point is too close to the epipole
        if ((ptMaxRound - ptEpipole2).squaredNorm() < 2500) continue;
        
        //commute the descriptor
        auto pt = round(toValidatePack.imagePointVec[idx]);
        CurveRasterizer<int, Polynomial2> descRaster(pt, epipoles.getFirstPx(),
                                            epipolarPtr->getFirst(toValidatePack.cloud[2*idx]));        
        if (epipoles.firstIsInverted()) descRaster.setStep(-1);
        vector<uint8_t> descriptor;
        const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
        
        Vector2i diff = round(ptMax - ptMin);
        const int distance = max(max(abs(diff[0]), abs(diff[1]))/step, 3) + MARGIN;
            
        CurveRasterizer<int, Polynomial2> raster(ptMaxRound, ptEpipole2,
                                            epipolarPtr->getSecond(toValidatePack.cloud[2*idx]));
        if (epipoles.secondIsInverted()) raster.setStep(-1);
        
        raster.setStep(step);
        raster.steps(-HALF_LENGTH);
        vector<uint8_t> sampleVec; //TODO declare outside the loop
        vector<int> uVec, vVec;
        uVec.reserve(distance);
        vVec.reserve(distance);
        sampleVec.reserve(distance);
        bool imageBorder = false;
        for (int d = 0; d < distance; d++, raster.step())
        {
            if (raster.v < 0 or raster.v >= img2.rows 
                or raster.u < 0 or raster.u >= img2.cols) 
                {
                    imageBorder = true;
                    break;
                }//sampleVec.push_back(0);
            else sampleVec.push_back(img2(raster.v, raster.u));
            uVec.push_back(raster.u);
            vVec.push_back(raster.v);
        }
        if (imageBorder) continue;
        
        //find the best correspondence
        vector<int> costVec = compareDescriptor(descriptor, sampleVec);
        auto bestCostIter = min_element(costVec.begin(), costVec.end());
        
        //check the cost
        Vector2i depthPt = depthValidPtVec[toValidatePack.idxMapVec[idx]];
        const int & x = depthPt[0], & y = depthPt[1];
        
        if (params.verbosity > 2)
        {
            cout << "Point : " << toValidatePack.imagePointVec[idx].transpose();
            cout << "   step : " << step << endl;
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
            cout << x << " " << y << "    " << *bestCostIter  << "    " <<  depth.cost(x, y) << endl << endl;
        }
        
        //TODO make sure that both stereo objects use the same descriptor length
        if (*bestCostIter > depth.cost(x, y) + 3*LENGTH or *bestCostIter > params.maxError)
        {
            //reject
            depth.at(x, y) = OUT_OF_RANGE;
            depth.sigma(x, y) = OUT_OF_RANGE;
        }
        else
        {
            //filter
            int dBest = bestCostIter - costVec.begin();
            double d1 = triangulate(toValidatePack.imagePointVec[idx][0], toValidatePack.imagePointVec[idx][1], 
                    uVec[dBest], vVec[dBest], CAMERA_1);
            double d2 = triangulate(toValidatePack.imagePointVec[idx][0], toValidatePack.imagePointVec[idx][1], 
                    uVec[dBest + 1], vVec[dBest + 1], CAMERA_1);
            double sigma1 = abs(d1 - d2) * SIGMA_COEFF;
            filter(depth.at(x, y), depth.sigma(x, y), d1, sigma1);
            depth.cost(x, y) = (depth.cost(x, y) + *bestCostIter)/2;
//            depth.at(x, y) = d1;
//            depth.sigma(x, y) = sigma1;
        }
    }
    
    
    ////for toCompute
    
    //simple reconstruction
    Vector3dVec toComputeCloud, toComputeCloud2;
    camera1->reconstructPointCloud(toComputePointVec, toComputeCloud);
    T12.inverseRotate(toComputeCloud, toComputeCloud2);
    cout << toComputePointVec.size() << endl;
    //search up to maximum disparity for the best correspondence
    for (int idx = 0; idx < toComputePointVec.size(); idx++)
    {
        
        // project min-max points
        Vector2d pt2;
        
        if (not camera2->projectPoint(toComputeCloud2[idx], pt2)) continue;
        
        const auto ptEpipole2 = epipoles.getSecondPx();
        const auto pt2Round = round(pt2);
        
        //the point is too close to the epipole
        if ((pt2Round - ptEpipole2).squaredNorm() < 2500) continue;
        
        const int step = stepVec[idx];
        const int distance = params.dispMax / step + MARGIN;
            
        CurveRasterizer<int, Polynomial2> raster(pt2Round, ptEpipole2,
                                            epipolarPtr->getSecond(toComputeCloud[idx]));
        if (epipoles.secondIsInverted()) raster.setStep(-1);
        
        raster.setStep(step);
        raster.steps(-HALF_LENGTH);
        vector<uint8_t> sampleVec; //TODO declare outside the loop
        vector<int> uVec, vVec;
        uVec.reserve(distance);
        vVec.reserve(distance);
        sampleVec.reserve(distance);
        bool imageBorder = false;
        for (int d = 0; d < distance; d++, raster.step())
        {
            if (raster.v < 0 or raster.v >= img2.rows 
                or raster.u < 0 or raster.u >= img2.cols) 
                {
                    imageBorder = true;
                    break;
                }//sampleVec.push_back(0);
            else sampleVec.push_back(img2(raster.v, raster.u));
            uVec.push_back(raster.u);
            vVec.push_back(raster.v);
        }
        if (imageBorder) continue;
        
        //find the best correspondence
        const vector<uint8_t> & descriptor = descriptorVec[idx]; 
        vector<int> costVec = compareDescriptor(descriptor, sampleVec);
        auto bestCostIter = min_element(costVec.begin(), costVec.end());
        
        //check the cost
        Vector2i depthPt = depthComputePtVec[idx];
        const int & x = depthPt[0], & y = depthPt[1];
        
        if (params.verbosity > 2)
        {
            cout << "Point : " << toComputePointVec[idx].transpose();
            cout << "   step : " << step << endl;
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
            cout << x << " " << y << "    " << *bestCostIter  << "    " <<  depth.cost(x, y) << endl << endl;
        }
        
        //TODO make sure that both stereo objects use the same descriptor length
        if (*bestCostIter < params.maxError)
        {
            //filter
            int dBest = bestCostIter - costVec.begin();
            double d1 = triangulate(toComputePointVec[idx][0], toComputePointVec[idx][1], 
                    uVec[dBest], vVec[dBest], CAMERA_1);
            double d2 = triangulate(toComputePointVec[idx][0], toComputePointVec[idx][1], 
                    uVec[dBest + 1], vVec[dBest + 1], CAMERA_1);
            double sigma1 = abs(d1 - d2) * SIGMA_COEFF;
            depth.at(x, y) = d1;
            depth.sigma(x, y) = sigma1;
            depth.cost(x, y) = *bestCostIter;
        }
    }
    
    //if the cost is OK, store the measurement
 
    delete epipolarPtr;     
}
