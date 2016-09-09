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

#pragma once


#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "utils/filter.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"

#include "reconstruction/stereo_misc.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/epipoles.h"
#include "reconstruction/eucm_stereo.h"

//TODO add errorMax threshold
struct MotionStereoParameters : StereoParameters
{
    int descLength = 5;
    int gradientThresh = 3;
};


class MotionStereo : private EnhancedStereo
{
public:
    MotionStereo(const EnhancedCamera * cam1, 
        const EnhancedCamera * cam2, MotionStereoParameters parameters) :
        EnhancedStereo(cam1, cam2, parameters),
        params(parameters)
    {
    }

    ~MotionStereo()
    {
        delete camera1;
        camera1 = NULL;
        delete camera2;
        camera2 = NULL;
    }    
    
    //TODO figure out how to treat the mask efficiently
    void setBaseImage(const Mat8u & image)
    {
        image.copyTo(img1);
        computeMask();
    }
       
    /*
    -Select salient points and points with defined depth
    -reproject all the points onto the next image
    -for those which have small uncertainties just keep that value
    -for those with wide uncertainty or without a value recompute it

    hypothesis quality must be evaluated using normal Kalman filtering
    That is, at every step the uncertainty grows. Once it reaches a certain value,
    the hypothesis is removed

    two sets of points must be treated separately:
    -points with bad descriptors but with depth estimation. 
        Project forward, increment uncertainty
    -points with good descriptors : 
        -with depth estimation. Project forward using serach 
        (if the search distance is greater than 1)
        -withoud depth estimation. Generate new hypotheses using complete epipolar search
    */
    //TODO cam1, cam2 are not needed
    void reprojectDepth(Transformation<double> T12, const Mat8u & img2, DepthMap & depth)
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
    
    void computeDepthOld(Transformation<double> T12,
            const Mat8u & img2, DepthMap & depth)
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
        EpipolarDescriptor epipolarDescriptor(LENGTH, WAVE_NORM, waveVec.data(), {1});
        
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
             QUERY_POINTS | DEFAULT_VALUES  | MINMAX | ALL_HYPOTHESES | SIGMA_VALUE | INDEX_MAPPING); 
        
        depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
        
        Vector3dVec cloud2;
        T12.inverseTransform(salientPack.cloud, cloud2);
        
        if (params.verbosity > 1) cout << "    searching for correspondences" << endl;
        for (int idx = 0; idx < salientPack.imagePointVec.size(); idx++)
        {
            
            // project min-max points
            Vector2d ptMin, ptMax;
            camera2->projectPoint(cloud2[2*idx], ptMin);
            camera2->projectPoint(cloud2[2*idx + 1], ptMax);
            
//            cout << cloud2[2*idx].norm() << " " << cloud2[2*idx + 1].norm() << endl;
            
            Vector2i diff = round(ptMax - ptMin);
            const int diffLength = max(abs(diff[0]), abs(diff[1])) + 1;
            if (diffLength > 3)
            {
                // if distance is big enough
                // search along epipolar curve
                
                // query 3D point must be projected into 1st frame
                CurveRasterizer<int, Polynomial2> raster(round(ptMax), round(ptMin),
                                                epipolarPtr->getSecond(salientPack.cloud[2*idx]));
                
                
                const int distance = min(diffLength, params.dispMax);
                
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
                
                //TODO make it possible to detect multiple hypotheses if there is no prior
                //TODO make this a parameter
                int dBest = -1;
                int eBest = LENGTH*15;
                for (int d = 0; d < distance; d++)
                {
                    const int & acc = costVec[d + HALF_LENGTH];
                    if (eBest > acc)
                    {
                        eBest = acc;
                        dBest = d;
                    }
                }
                if (dBest == -1) continue;
                // triangulate and improve sigma
                double d1 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                        uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], CAMERA_1);
                double d2 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                        uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], CAMERA_1);
                double sigma1 = abs(d2 - d1) / 1.732; //variance of a uniform distribution
                
                //TODO make push by idx
                Vector2i depthPt = depthPointVec[salientPack.idxMapVec[idx]];
                depth.pushHypothesis(depthPt[0], depthPt[1], d1, sigma1);
            }
            
        }
    }
    
    void computeDepth(Transformation<double> T12,
            const Mat8u & img2, DepthMap & depth)
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
        EpipolarDescriptor epipolarDescriptor(LENGTH, WAVE_NORM, waveVec.data(), {1});
        
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
             QUERY_POINTS | MINMAX | ALL_HYPOTHESES | SIGMA_VALUE | INDEX_MAPPING); 
        
        depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
        
        Vector3dVec cloud2;
        T12.inverseTransform(salientPack.cloud, cloud2);
        
        if (params.verbosity > 1) cout << "    searching for correspondences" << endl;
        for (int idx = 0; idx < salientPack.imagePointVec.size(); idx++)
        {
            
            // project min-max points
            Vector2d ptMin, ptMax;
            camera2->projectPoint(cloud2[2*idx], ptMin);
            camera2->projectPoint(cloud2[2*idx + 1], ptMax);
            
//            cout << cloud2[2*idx].norm() << " " << cloud2[2*idx + 1].norm() << endl;
            
            Vector2i diff = round(ptMax - ptMin);
            const int diffLength = max(abs(diff[0]), abs(diff[1])) + 1;
            if (diffLength > 3)
            {
                // if distance is big enough
                // search along epipolar curve
                
                // query 3D point must be projected into 1st frame
                CurveRasterizer<int, Polynomial2> raster(round(ptMax), round(ptMin),
                                                epipolarPtr->getSecond(salientPack.cloud[2*idx]));
                
                
                const int distance = diffLength;
                
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
                
                //TODO make it possible to detect multiple hypotheses if there is no prior
                //TODO make this a parameter
                int dBest = -1;
                int eBest = LENGTH*15;
                for (int d = 0; d < distance; d++)
                {
                    const int & acc = costVec[d + HALF_LENGTH];
                    if (eBest > acc)
                    {
                        eBest = acc;
                        dBest = d;
                    }
                }
                if (dBest == -1) continue;
                // triangulate and improve sigma
                double d1 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                        uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], CAMERA_1);
                double d2 = triangulate(salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], 
                        uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], CAMERA_1);
                double sigma1 = abs(d2 - d1) / 1.732; //variance of a uniform distribution
                
                //TODO make push by idx
                Vector2i depthPt = depthPointVec[salientPack.idxMapVec[idx]];
                
                depth.at(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = d1;
                depth.sigma(depthPt[0], depthPt[1], salientPack.hypIdxVec[idx]) = sigma1;
            }
            
        }
        depth.regularize();
    }
    
private:
    
    EnhancedEpipolar * epipolarPtr;
    // based on the image gradient
    void computeMask()
    {
        Mat16s gradx, grady;
        Sobel(img1, gradx, CV_16S, 1, 0, 1);
        Sobel(img1, grady, CV_16S, 0, 1, 1);
        Mat16s gradAbs = abs(gradx) + abs(grady);
        GaussianBlur(gradAbs, gradAbs, Size(7, 7), 0, 0);
        Mat8u gradAbs8u;
        gradAbs.convertTo(gradAbs8u, CV_8U);
        threshold(gradAbs8u, maskMat, params.gradientThresh, 128, CV_THRESH_BINARY);
        
    }
   
    Mat8u img1;    
    Mat8u maskMat; //TODO compute mask
    const MotionStereoParameters params;
};

