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


//TODO add errorMax threshold
struct MotionStereoParameters
{
    int scale = 3;
    int descLength = 5;
    int gradientThresh = 3;
    int verbosity = 0;
    int uMargin = 25, vMargin = 25;  // RoI left upper corner
    int biasMax = 10;
    int dispMax = 25;
};


class MotionStereo
{
public:
    MotionStereo(const EnhancedCamera * cam1, 
        const EnhancedCamera * cam2, MotionStereoParameters params) :
        // initialize the members
        camera1(cam1->clone()),
        camera2(cam2->clone()),
        params(params)
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
    
    //FIXME temporary function
    const int FLAW_COST = 7;
    vector<int> compareDescriptor(const vector<uint8_t> & desc, const vector<uint8_t> & sampleVec)
    {
        const int HALF_LENGTH = desc.size() / 2;
        vector<int> rowA(sampleVec.size()), rowB(sampleVec.size());
        
        //match the first half
        for (int i = 0; i < sampleVec.size(); i++)
        {
            rowA[i] = (abs(sampleVec[i] - desc[0]));
        }
        for (int i = 1; i <= HALF_LENGTH; i++)
        {
            rowB[0] = (rowA[0] + FLAW_COST + abs(sampleVec[0] - desc[i]));
            int cost = min(rowA[i] + FLAW_COST, rowA[0]);
            rowB[1] = (cost + abs(sampleVec[i] - desc[i]));
            for (int j = 2; j < sampleVec.size(); j++)
            {
                cost = min(min(rowA[j] + FLAW_COST, rowA[j - 1]), rowA[j - 2] + FLAW_COST);
                rowB[j] = (cost+ abs(sampleVec[j] - desc[i]));
            }
            swap(rowA, rowB);
        }
        vector<int> rowC(sampleVec.size()); //center cost
        swap(rowA, rowC);
        
        //match the second half (from the last pixel to first)
        for (int i = 0; i < sampleVec.size(); i++)
        {
            rowA[i] = (abs(sampleVec[i] - desc.back()));
        }
        for (int i = desc.size() - 1; i > HALF_LENGTH + 1; i--)
        {
            for (int j = 0; j < sampleVec.size() - 2; j++)
            {
                int cost = min(min(rowA[j] + FLAW_COST, rowA[j + 1]), rowA[j + 2] + FLAW_COST);
                rowB[j] = (cost + abs(sampleVec[j] - desc[i]));
            }
            int j = sampleVec.size() - 2;
            int cost = min(rowA[j] + FLAW_COST, rowA[j + 1]);
            rowB[j] = (cost + abs(sampleVec[j] - desc[i]));
            rowB.back() = (rowA.back() + FLAW_COST + abs(sampleVec.back() - desc[i]));
            swap(rowA, rowB);
        }
        
        //accumulate the cost
        for (int i = 0; i < sampleVec.size() - 2; i++)
        {
            rowC[i] += min(min(rowA[i] + FLAW_COST, rowA[i + 1]), rowA[i + 2] + FLAW_COST);
        }
        int i = rowC.size() - 2;
        rowC[i] += min(rowA[i] + FLAW_COST, rowA[i + 1]);
        rowC.back() += rowA.back() + FLAW_COST;
        return rowC;
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
                    //FIXME if some points are not reconstructed indexes are not matched
                    descriptorVec.push_back(descriptor);
                    salientPack.imagePointVec.push_back(pt);
                }
            }
        }
        
        //TODO check that sigmaVec is reconstructed
        depth.reconstruct(flatPack, QUERY_POINTS | SIGMA_VALUE);
        depth.reconstruct(salientPack, QUERY_POINTS /*| DEFAULT_VALUES*/  | MINMAX | ALL_HYPOTHESES | SIGMA_VALUE | INDEX_MAPPING);
        depth.setTo(OUT_OF_RANGE, OUT_OF_RANGE);
        
        // for the flat pack project points and replace the hypotheses in depth
        
        T12.inverseTransform(flatPack.cloud, flatPack.cloud);
//        cout << flatPack.cloud.size() << "  " << flatPack.sigmaVec.size() << endl;
        for (int idx = 0; idx < flatPack.cloud.size(); idx++)
        {
            depth.pushHypothesis(flatPack.cloud[idx], flatPack.sigmaVec[idx]);
        }
        
        // for the salient pack compute stereo and for corresponding pixel push new hypothesis
        Vector3dVec cloud2;
        T12.inverseTransform(salientPack.cloud, cloud2);
        Transform12 = T12.inverse();
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
            else
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
                
                int dBest = 0;
                int eBest = LENGTH*65535;
                int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
                for (int d = 0; d < distance; d++)
                {
                    int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
                    int bias = 0; //min(params.biasMax, max(-params.biasMax, (sum2 - sum1) / NORMALIZER));
                    int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
                                        descriptor.begin(), sampleVec.begin() + d, bias, 1);
                    if (eBest > acc)
                    {
                        eBest = acc;
                        dBest = d;
                    }
                }
                dispAcc += dBest;
                dispCount++;
//                auto poly1 = epipolarPtr->getSecond(salientPack.cloud[2*idx]);
//                
//                cout << "    curve : " << poly1(ptMin[0], ptMin[1]) << " " 
//                    << poly1(ptMax[0], ptMax[1]) <<  " " << poly1(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH]) <<endl;
                
//                cout << ptMax.transpose() << " / " << uVec[dBest + HALF_LENGTH] << " " 
//                    << vVec[dBest + HALF_LENGTH] << " / " << ptMin.transpose() << " / " << dBest << " " << distance << endl;
                
                // put the original depth to the new pixel
           
//                int xd = depth.xConv(uVec[dBest + HALF_LENGTH]);
//                int yd = depth.yConv(vVec[dBest + HALF_LENGTH]);
//                if (not depth.isValid(xd, yd)) continue;
//                depth.at(xd, yd) = 0.5*(salientPack.cloud[2*idx] + salientPack.cloud[2*idx + 1]).norm();
//                depth.sigma(xd, yd) = 1;

                // put the original hypothesis
//                depth.pushHypothesis(0.5*(salientPack.cloud[2*idx] + salientPack.cloud[2*idx + 1]), 1);

                // triangulate and improve sigma
                 Vector3d X1, X2;
                triangulate(uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], 
                        salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], X1);
                triangulate(uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], 
                        salientPack.imagePointVec[idx][0], salientPack.imagePointVec[idx][1], X2);
//                cout << X1.transpose() - 0.5*(salientPack.cloud[2*idx] + salientPack.cloud[2*idx + 1]).transpose() << endl;
                depth.pushHypothesis(X1, (X2 - X1).norm() / 2);
            }
            
        }     
        cout << double(dispAcc) / dispCount << endl;
        //release dynamic objects
        delete epipolarPtr;
        epipolarPtr = NULL;
    }    
    
    //TODO split into functions
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
        EpipolarDescriptor epipolarDescriptor(LENGTH, WAVE_NORM, waveVec.data(), {1, 2, 3});
        
        if (params.verbosity > 1) cout << "    computing the scan limits" << endl;
        // get uncertainty range reconstruction in the first frame
        
        // discard non-salient points
        depth.applyMask(maskMat);
        
        vector<int> idxVec;
        Vector3dVec minDistVec, maxDistVec;
        depth.reconstructUncertainty(idxVec, minDistVec, maxDistVec);
        Vector2dVec pointVec = depth.getPointVec(idxVec);
        // reproject them onto the second image
        Vector3dVec minDist2Vec, maxDist2Vec;
        T12.inverseTransform(minDistVec, minDist2Vec);
        T12.inverseTransform(maxDistVec, maxDist2Vec);
        Vector2dVec pointMinVec;
        Vector2dVec pointMaxVec;
        vector<bool> maskMinVec, maskMaxVec;
        camera2->projectPointCloud(minDist2Vec, pointMinVec, maskMinVec);
        camera2->projectPointCloud(maxDist2Vec, pointMaxVec, maskMaxVec);
        
        vector<bool> maskVec;
        for (int i = 0; i < maskMinVec.size(); i++)
        {
            maskVec.push_back(maskMinVec[i] and maskMaxVec[i]);
        }
        
        
        if (params.verbosity > 1) cout << "    core loop" << endl;
        for (int ptIdx = 0; ptIdx < minDistVec.size(); ptIdx++)
        {
            if (not maskVec[ptIdx])
            {
                depth.at(idxVec[ptIdx]) = 0;
                continue;
            }   
            
            // ### compute descriptor ###
            if (params.verbosity > 2) cout << "        compute descriptor" << endl;
            // get the corresponding rasterizer
            CurveRasterizer<int, Polynomial2> descRaster(round(pointVec[ptIdx]), epipoles.getFirstPx(),
                                                epipolarPtr->getFirst(minDistVec[ptIdx]));
            if (epipoles.firstIsInverted()) descRaster.setStep(-1);
            
            // compute the actual descriptor
            vector<uint8_t> descriptor;
            const int step = epipolarDescriptor.compute(img1, descRaster, descriptor);
            if (not epipolarDescriptor.goodResp() or step < 1)
            {
                depth.at(idxVec[ptIdx]) = 0;
                continue;
            }
            // ### find the best correspondence on img2 ###
            if (params.verbosity > 2) cout << "        sampling the second image" << endl;
            //sample the second image
            //TODO traverse the epipolar line in the opposit direction and respect the disparity limit
            Vector2i goal = round(pointMinVec[ptIdx]);
            Vector2i start = round(pointMaxVec[ptIdx]);
            CurveRasterizer<int, Polynomial2> raster(start, goal,
                                                epipolarPtr->getSecond(minDistVec[ptIdx]));
            Vector2i diff = start - goal;
            
            /*{
                Vector3d X1, X2;
                triangulate(pointVec[ptIdx][0], pointVec[ptIdx][1], 
//                        pointMinVec[ptIdx][0], pointMinVec[ptIdx][1], X1);
                        goal[0], goal[1], X1);
                triangulate(pointVec[ptIdx][0], pointVec[ptIdx][1], 
//                        pointMaxVec[ptIdx][0], pointMaxVec[ptIdx][1], X2);
                        start[0], start[1], X2);
                cout << "MIN MAX" << endl;
                cout << minDistVec[ptIdx].transpose() - X1.transpose() << endl;
                cout << maxDistVec[ptIdx].transpose() - X2.transpose() << endl;
            }*/
            
            
            const int distance = min(int(diff.norm()), params.dispMax);
            raster.steps(-HALF_LENGTH*step);
            vector<uint8_t> sampleVec;
            vector<int> uVec, vVec;
            const int margin = LENGTH - 1;
            uVec.reserve(distance + margin);
            vVec.reserve(distance + margin);
            sampleVec.reserve(distance + margin);
            for (int d = 0; d < distance + margin; d++, raster.steps(step))
            {
                if (raster.v < 0 or raster.v >= img2.rows 
                    or raster.u < 0 or raster.u >= img2.cols) sampleVec.push_back(0);
                else sampleVec.push_back(img2(raster.v, raster.u));
                uVec.push_back(raster.u);
                vVec.push_back(raster.v);
            }
            
            if (params.verbosity > 2) cout << "        find the best candidate" << endl;
            
            vector<int> costVec = compareDescriptor(descriptor, sampleVec);
            
            //compute the error and find the best candidate
            int dBest = 0;
            int eBest = LENGTH*255;
//            int sum1 = filter(kernelVec.begin(), kernelVec.end(), descriptor.begin(), 0);
//            cout << "ERROR CURVE " << step << endl;
            for (int d = 0; d < distance; d++)
            {
//                int sum2 = filter(kernelVec.begin(), kernelVec.end(), sampleVec.begin() + d, 0);
//                int bias = min(params.biasMax, max(-params.biasMax, (sum2 - sum1) / NORMALIZER));
//                int acc =  biasedAbsDiff(kernelVec.begin(), kernelVec.end(),
//                                    descriptor.begin(), sampleVec.begin() + d, bias, step);
//                cout << acc << endl;
                int acc = costVec[d + HALF_LENGTH];
                if (eBest > acc)
                {
                    eBest = acc;
                    dBest = d;
                }
            }
            
            //TODO make triangulation checks and 
            
            Vector3d X1, X2;
            triangulate(pointVec[ptIdx][0], pointVec[ptIdx][1], 
                    uVec[dBest + HALF_LENGTH], vVec[dBest + HALF_LENGTH], X1);
            triangulate(pointVec[ptIdx][0], pointVec[ptIdx][1], 
                    uVec[dBest + HALF_LENGTH + 1], vVec[dBest + HALF_LENGTH + 1], X2);
            depth.at(idxVec[ptIdx]) = X1.norm();
            depth.sigma(idxVec[ptIdx]) = (X2 - X1).norm()/2;
//            if (distance > 5)
//            {
//                cout << "BAD DISTANCE: " << distance << endl;
//                cout << start.transpose() << "   " << goal.transpose() << endl;        
//            }
//            if ( depth.sigma(idxVec[ptIdx]) > 0.5)
//            {
//                cout << distance << endl;
//                cout << "depth: " << depth.at(idxVec[ptIdx]) << " +-" << depth.sigma(idxVec[ptIdx]) << endl;
//                cout << "samples:" << endl;
//                for (auto & x : sampleVec)
//                {
//                    cout << " " << int(x);
//                }
//                cout << endl;
//                cout << "descriptor:" << endl;
//                for (auto & x : descriptor)
//                {
//                    cout << " " << int(x);
//                }
//                cout << endl;
//                cout << "cost:" << endl;
//                for (auto & x : costVec)
//                {
//                    cout << " " << int(x);
//                }
//                cout << endl<< endl;
//            }
            
//            cout << minDistVec[ptIdx].norm() << "  " << X1.norm() 
//                    << "  " << maxDistVec[ptIdx].norm() << endl;
        }
        
        delete epipolarPtr;
        epipolarPtr = NULL;
    }
    
   
    // TODO a lot of overlap with EnhancedStereo, think about merging them or deriving them
    bool triangulate(double x1, double y1, double x2, double y2, Vector3d & X)
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
    
    // pose of the first to the second camera
    Transformation<double> Transform12;
    EnhancedCamera *camera1, *camera2;
    
    Mat8u img1;    
    Mat8u maskMat; //TODO compute mask
    const MotionStereoParameters params;
};

