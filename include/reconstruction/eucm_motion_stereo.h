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
#include "filter.h"
#include "geometry/geometry.h"
#include "camera/eucm.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"

struct MotionStereoParameters
{
    int descLength = 5;
    int gradientThersh = 5;
    int verbosity = 0;
    int uMargin = 25, vMargin = 25;  // RoI left upper corner
    int maxBias = 10;
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
    
    void setBaseImage(const Mat8u & image);
    
    void computeDepth(Transformation<double> T12, const DepthMap & depthPrior,
            Mat8u img2, DepthMap & depthOut)
    {
        epipolarPtr = new EnhancedEpipolar(T12, camera1, camera2, 2000, verbosity);
        Transform12 = T12;
        // init the output mat
        //TODO think about overhead
        depthOut = depthPrior;
        depthOut.setTo(0, 1);
        
        int LENGTH = params.descLength;
        int HALF_LENGTH = LENGTH / 2;
        
        // compute the weights for matching cost
        // TODO make a separate header and a cpp file miscellaneous
        vector<int> kernelVec(LENGTH);
        int NORMALIZER;
        switch (LENGTH)
        {
        case 3:
            copy(KERNEL_3.begin(), KERNEL_3.end(), kernelVec.begin());
            NORMALIZER = NORMALIZER_3;
            break;
        case 5:
            copy(KERNEL_5.begin(), KERNEL_5.end(), kernelVec.begin());
            NORMALIZER = NORMALIZER_5;
            break;
        case 7:
            copy(KERNEL_7.begin(), KERNEL_7.end(), kernelVec.begin());
            NORMALIZER = NORMALIZER_7;
            break;
        default:
            LENGTH = 9;
            HALF_LENGTH = 4;
        case 9:
            copy(KERNEL_9.begin(), KERNEL_9.end(), kernelVec.begin());
            NORMALIZER = NORMALIZER_9;
            break;
        }
        
        // get uncertainty range reconstruction in the first frame
        Vector2dVec pointVec; 
        Vector3dVec minDistVec, maxDistVec;
        depthPrior.reconstructUncertainty(pointVec, minDistVec, maxDistVec);
        
        // discard non-salient points
        vector<bool> maskVec;
        for (auto & pt : pointVec)
        {
            if ( maskMat(round(pt[1]), round(pt[0])) ) maskVec.push_back(true);
            else maskVec.push_back(false);
        }
        
        // reproject them onto the second image
        Vector3dVec minDist2Vec, maxDist2Vec;
        T12.inverseTransform(minDistVec, minDist2Vec);
        T12.inverseTransform(maxDistVec, maxDist2Vec);
        
        Vector2dVec pointMinVec;
        Vector2dVec pointMaxVec;
        vector<bool> maskMinVec, maskMaxVec;
        camera2->projectPointCloud(minDist2Vec, pointMinVec, maskMinVec);
        camera2->projectPointCloud(maxDist2Vec, pointMaxVec, maskMaxVec);
        
        for (int i = 0; i < maskVec.size(); i++)
        {
            if (not maskMinVec[i] or not maskMaxVec[i]) maskVec[i] = false;
        }
        
        
        for (int ptIdx = 0; ptIdx < minDistVec.size(); ptIdx++)
        {
            if (not maskVec[ptIdx])
            {
                depthOut.nearest(pointVec[ptIdx]) = 0;
                continue;
            }   
            
            // ### compute descriptor ###
            
            // get the corresponding rasterizer
            CurveRasterizer<int, Polynomial2> descRaster(round(pointVec[ptIdx]), epipolePx1,
                                                epipolarPtr->getFirst(minDistVec[ptIdx]));
            if (epipoleInverted1) descRaster.setStep(-1);
            
            // compute the actual descriptor
            descRaster.steps(-HALF_LENGTH);
            vector<uint8_t> descriptor(LENGTH);
            for (int i = 0; i < LENGTH; i++, descRaster.step())
            {
                descriptor[i] = img1(descRaster.v, descRaster.u);
            }
            
            // ### find the best correspondence on img2 ###
            
            //sample the second image
            CurveRasterizer<int, Polynomial2> raster(round(pointMinVec[ptIdx]), round(pointMaxVec[ptIdx]),
                                                epipolarPtr->getSecond(minDistVec[ptIdx]));
            raster.steps(-HALF_LENGTH);
            vector<uint8_t> sampleVec;
            Vector2i goal = round(pointMaxVec[ptIdx]);
            vector<int> uVec, vVec;
            int loopCount = HALF_LENGTH;
            bool loopFlag = false;
            while (loopCount)
            {
                if (raster.v < 0 or raster.v >= img2.rows 
                    or raster.u < 0 or raster.u >= img2.cols) sampleVec.push_back(0);
                else sampleVec.push_back(img2(raster.v, raster.u));
                uVec.push_back(raster.u);
                vVec.push_back(raster.v);
                raster.step();
                if (not loopFlag) loopFlag = (abs(goal[0] - raster.u) + abs(goal[1] - raster.v) <= 1);
                else loopCount--;
            }
            
            //compute the error and find the best candidate
            int dBest = 0;
            int eBest = LENGTH*255;
            int sum1 = accumulate(descriptor.begin(), descriptor.end(), 0);
            for (int d = 0; d < sampleVec.size() - LENGTH + 1; d++)
            {
                int acc = 0;
                int sum2 = accumulate(sampleVec.begin() + d, sampleVec.begin() + d + LENGTH, 0);
                int bias = min(params.maxBias, max(-params.maxBias, (sum2 - sum1) / LENGTH));
                for (int i = 0; i < LENGTH; i++)
                {
                    acc += abs(descriptor[i] - sampleVec[d + i] + bias) * kernelVec[i];
                }
                
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
            depthOut.nearest(pointVec[ptIdx]) = X1.norm();
            depthOut.nearestSigma(pointVec[ptIdx]) = (X2.norm() - X1.norm()) / 2;
        }
        
        delete epipolarPtr;
        epipolarPtr = NULL;
    }
    
   
    // TODO a lot of overlap with EnhancedStereo, think about merging them of deriving them
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
    
    void computeEpipole(Transformation<double> Transform12)
    {
        if (not camera1->projectPoint(Transform12.trans(), epipole1))
        {
            camera1->projectPoint(-Transform12.trans(), epipole1);
            epipoleInverted1 = true;
        }
        else epipoleInverted1 = false;
        
        //FIXME not nedeed?
        if (not camera2->projectPoint(Transform12.transInv(), epipole2))
        {
            camera2->projectPoint(-Transform12.transInv(), epipole2);
            epipoleInverted2 = false;
        }
        else epipoleInverted2 = false;
        
        epipolePx1 = round(epipole1);
        epipolePx2 = round(epipole2);
    }
    
private:
    EnhancedEpipolar * epipolarPtr;
    // based on the image gradient
    void computeMask();
    
    bool epipoleInverted1, epipoleInverted2;
    Vector2d epipole1, epipole2;  // projection of the first camera center onto the second camera
    Vector2i epipolePx1, epipolePx2;
    
    // pose of the first to the second camera
    Transformation<double> Transform12;
    EnhancedCamera *camera1, *camera2;
    
    Mat8u img1;    
    Mat8u maskMat; //TODO compute mask
    MotionStereoParameters params;
    int verbosity;
};

