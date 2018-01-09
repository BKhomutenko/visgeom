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
Relative camera pose estimation based on photometric error and depth map
*/

#include "localization/sparse_odom.h"

#include "std.h"
#include <tuple>
#include "eigen.h"
#include "ocv.h"
#include "ceres.h" 

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "reconstruction/triangulator.h"
#include "localization/local_cost_functions.h"

using std::get;
using std::tie;
using std::tuple;

using cv::Vec3b;
int imgCount = 0;



void harrisResp(const Mat8u & img, Mat32f & respMat)
{
    Mat32f gx, gy;
    Mat32f Ixx, Ixy, Iyy;
    Sobel(img, gx, CV_32F, 1, 0, 3, 0.01);
    Sobel(img, gy, CV_32F, 0, 1, 3, 0.01);
    Ixx = gx.mul(gx);
    Ixy = gx.mul(gy);
    Iyy = gy.mul(gy);
    GaussianBlur(Ixx, Ixx, Size(3, 3), 0.5, 0.5);
    GaussianBlur(Ixy, Ixy, Size(3, 3), 0.5, 0.5);
    GaussianBlur(Iyy, Iyy, Size(3, 3), 0.5, 0.5);
    respMat = (Ixx.mul(Iyy) - Ixy.mul(Ixy) )/ ( 1 + Ixx + Iyy);
    GaussianBlur(respMat, respMat, Size(3, 3), 0.7, 0.7);
}

void harrisBlobs(const Mat32f & respMat, vector<pair<double, int>> respHeap, const Vector2dVec &  maxVec)
{
    
    Mat8u visited(respMat.size());
    
    Mat8uc3 visited2(respMat.size());
    visited2.setTo(0);
    imshow("harris", respMat);
    int countFeatures = 0;
    while (countFeatures < 100 and not respHeap.empty())
    {
        pop_heap(respHeap.begin(), respHeap.end());
        
        vector<tuple<double, int, int>> heap;
        
        int idx = respHeap.back().second;
        double resp = respHeap.back().first;
        respHeap.pop_back();
        if (visited2(maxVec[idx][1], maxVec[idx][0])[0] != 0) continue;
        countFeatures++;
        heap.emplace_back(respMat(maxVec[idx][1], maxVec[idx][0]), maxVec[idx][0], maxVec[idx][1]);
        visited.setTo(0);
        double acc = 0;
        int count = 0;
        while (not heap.empty())
        {
            pop_heap(heap.begin(), heap.end());
            double resp;
            int u, v;
            tie(resp, u, v) = heap.back();
            heap.pop_back();
            if (visited2(v, u)[0] != 0 or visited(v, u) == 255) continue;
            if (2 * count * resp < acc) continue;
            visited(v, u) = 255;
            visited2(v, u) = Vec3b((countFeatures*17) % 255 + 150,
                                    (countFeatures*37) % 255 + 150,
                                    (countFeatures*47) % 255 + 150);
                                    
            acc += resp;
            count++;
            if (u > 0 and 2 * count * respMat(v, u - 1) > acc)
            {
                heap.emplace_back(respMat(v, u - 1), u - 1, v);
                push_heap(heap.begin(), heap.end());
            }
            if (u < respMat.cols - 1 and 2 * count * respMat(v, u + 1) > acc)
            {
                heap.emplace_back(respMat(v, u + 1), u + 1, v);
                push_heap(heap.begin(), heap.end());
            }
            if (v > 0 and 2 * count * respMat(v - 1, u) > acc)
            {
                heap.emplace_back(respMat(v - 1, u), u, v - 1);
                push_heap(heap.begin(), heap.end());
            }
            if (v < respMat.rows - 1 and 2 * count * respMat(v + 1, u) > acc)
            {
                heap.emplace_back(respMat(v + 1, u), u, v);
                push_heap(heap.begin(), heap.end());
            }
        }
//        
    }  
    imshow("visited" + to_string(imgCount++), visited2);
    waitKey();
}

void harrisTest()
{
    Mat32f respMat;
    Mat8u img = imread("/home/bogdan/projects/stack/lena.png", 0);
    harrisResp(img, respMat);
    vector<pair<double, int>> respHeap;
    Vector2dVec maxVec;
    
    for (int v = 7; v < img.rows - 7; v++)
    {
        for (int u = 7; u < img.cols - 7; u++)
        {
            double respVal = respMat(v, u);
            bool isMax = true;
            for (int dv = -1; dv <= 1 and isMax; dv++)
            {
                for (int du = -1; du <= 1 and isMax; du++)
                {
                    if (du == 0 and dv == 0) continue;
                    if (respVal <= respMat(v + dv, u + du)) isMax = false;
                }
            }
            if (isMax)
            {
                maxVec.emplace_back(u, v);
                respHeap.emplace_back(respVal, maxVec.size() - 1);
            }
        }
    }
    make_heap(respHeap.begin(), respHeap.end());
    harrisBlobs(respMat, respHeap, maxVec);
}

//FIXME TMP HARRIS FEATURES
Vector2dVec harrisCorners(const Mat8u & img)
{
//    harrisTest();
    Mat32f resp;
    resp.create(img.size());
    cornerHarris(img, resp, 7, 3, 0.05);
//    harrisResp(img, resp);
    vector<pair<double, int>> respHeap;
    Vector2dVec maxVec;
    
    for (int v = 7; v < img.rows - 7; v++)
    {
        for (int u = 7; u < img.cols - 7; u++)
        {
            double respVal = resp(v, u);
            bool isMax = true;
            for (int dv = -1; dv <= 1 and isMax; dv++)
            {
                for (int du = -1; du <= 1 and isMax; du++)
                {
                    if (du == 0 and dv == 0) continue;
                    if (respVal <= resp(v + dv, u + du)) isMax = false;
                }
            }
            if (isMax)
            {
                maxVec.emplace_back(u, v);
                respHeap.emplace_back(respVal, maxVec.size() - 1);
            }
        }
    }
    make_heap(respHeap.begin(), respHeap.end());
    Vector2dVec resVec;
    const int NUM_FEATURES = 500;
    while (resVec.size() < NUM_FEATURES and not respHeap.empty())
    {
        pop_heap(respHeap.begin(), respHeap.end());
        int idx = respHeap.back().second;
        respHeap.pop_back();
        resVec.push_back(maxVec[idx]);
    } 
    return resVec;
}

void descriptors(const Mat8u & img, Mat32f & out, const Vector2dVec & pointVec, const int DESC_SIZE = 4)
{
    out.create( Size(pow(2 * DESC_SIZE + 1, 2), pointVec.size()) );
    vector<pair<double, int>> respHeap;
    Vector2dVec maxVec;
    
    vector<double> kernel;
    
    for (int v = -DESC_SIZE; v <= DESC_SIZE; v++)
    {
        double y = (2. * v) / DESC_SIZE;
        for (int u = -DESC_SIZE; u <= DESC_SIZE; u++)
        {   
            double x = (2. * u) / DESC_SIZE;
            kernel.push_back( exp(-0.5*(x*x + y*y)) );
        }
    }
    
    for (int idx = 0; idx < pointVec.size(); idx++)
    {
        bool isMax = true;
        const int u = round(pointVec[idx][0]);
        const int v = round(pointVec[idx][1]);
        float * outPtr = (float*)(out.row(idx).data);
        auto kernelIter = kernel.begin();
        for (int dv = -DESC_SIZE; dv <= DESC_SIZE; dv++)
        {
            for (int du = -DESC_SIZE; du <= DESC_SIZE; du++, outPtr++, kernelIter++)
            {   
                *outPtr = *kernelIter * img(v + dv, u + du);
            }
        }
    }
}

void SparseOdometry::feedData(const Mat8u & imageNew, const Transf xiOdomNew)
{
    Transf dxi = xiBaseCam.inverseCompose(xiOdom.inverseCompose(xiOdomNew)).compose(xiBaseCam);
    //FIXME for debug
    if (dxi.trans().norm() < MIN_STEREO_BASE and not keypointVec1.empty()) return;
    
    cout << "DETECT" << endl;
//    keypointVec2 = harrisCorners(imageNew.rowRange(0, 300));
    keypointVec2 = harrisCorners(imageNew);
    cout << "EXTRACT" << endl;
    descriptors(imageNew, desc2, keypointVec2);
    
//    detector.detect(imageNew.rowRange(0, 300), keypointVec2);
    
//    for (auto & x : keypointVec2)
//    {
//        x.angle = 0;
//    }
    
//    detector.compute(imageNew, keypointVec2, desc2);

//    cout << "detected : " << keypointVec2.size() << endl; 
    if (not keypointVec1.empty())
    {
        cout << "MATCH" << endl;
        //correspondence
        BFMatcher matcher(cv::NORM_L1, true);
        vector<DMatch> matchVec;
        matcher.match(desc1, desc2, matchVec);
        
        
        
        
        
        
        
        
        //reconstruct
        vector<DMatch> goodMatchVec;
        Vector2dVec pointVec1, pointVec2;
        vector<double> sizeVec;
        for (auto & match : matchVec)
        {
            if (match.distance > distThresh) continue;
//            if (abs(keypointVec1[match.queryIdx].size - keypointVec2[match.trainIdx].size) /
//                (keypointVec1[match.queryIdx].size + keypointVec2[match.trainIdx].size) > 0.07) continue;
            goodMatchVec.push_back(match);
            
//            const auto & pt1 = keypointVec1[match.queryIdx].pt;
//            const auto & pt2 = keypointVec2[match.trainIdx].pt;
//            
////                if (norm(pt1 - pt2) < 3) continue;
//            
//            pointVec1.emplace_back(pt1.x, pt1.y);
//            pointVec2.emplace_back(pt2.x, pt2.y);

            pointVec1.push_back(keypointVec1[match.queryIdx]);
            pointVec2.push_back(keypointVec2[match.trainIdx]);
            
            sizeVec.push_back(1);
        }
//        cout << endl;
//        Mat img_matches;
//        vector<KeyPoint> kp1Vec, kp2Vec;
//        for (auto & x : keypointVec1)
//        {
//            kp1Vec.emplace_back(x[0], x[1], 1);
//        }
//        for (auto & x : keypointVec2)
//        {
//            kp2Vec.emplace_back(x[0], x[1], 1);
//        }
//        drawMatches( imageOld, kp1Vec, imageNew, kp2Vec,
//            goodMatchVec, img_matches, Scalar::all(-1), Scalar::all(-1),
//            vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
//        imshow("matches", img_matches);
//        waitKey();
//        
        

        Vector3dVec reconstVec1, reconstVec2;
        camera->reconstructPointCloud(pointVec1, reconstVec1);
        camera->reconstructPointCloud(pointVec2, reconstVec2);
        
        
        
        
        
        //RANSAC
        cout << "RANSAC" << endl;
        vector<bool> inlierMask;
        ransacNPoints(reconstVec1, reconstVec2, pointVec2, sizeVec,
                xiOdom.inverseCompose(xiOdomNew), inlierMask);
        
        
        
        cout << "FINAL" << endl;
        //refinement
        Vector3dVec xInlierVec1, xInlierVec2;
        Vector2dVec pInlierVec2;
        vector<double> sizeVec2;
        for (int i = 0; i < inlierMask.size(); i++)
        {
            if (not inlierMask[i]) continue;
            
            xInlierVec1.push_back(reconstVec1[i]);
            xInlierVec2.push_back(reconstVec2[i]);
            pInlierVec2.push_back(pointVec2[i]);
            sizeVec2.push_back(sizeVec[i]);
        }
        computeTransfSparse(xInlierVec1, xInlierVec2, pInlierVec2, sizeVec2,
                xiOdom.inverseCompose(xiOdomNew), xiIncr, true);
                
        ///// Recompute the motions with discarded close outliers
        Transf xiTest = xiBaseCam.inverseCompose(xiIncr).compose(xiBaseCam);
        Vector2dVec reprojectedVec = reprojectPoints(xInlierVec1, xInlierVec2, dxi);
        
        // compute sigma
        double sigmaSqAcc = 0;
        vector<double> errVec;
        for (int i = 0; i < reprojectedVec.size(); i++)
        {
            const double err = (pInlierVec2[i] - reprojectedVec[i]).squaredNorm();
            errVec.push_back(err);
            sigmaSqAcc += err;
        }
        const double sigmaSq = sigmaSqAcc / (reprojectedVec.size() - 2); // 2 degrees of freedom
        
        Vector3dVec xInlierVec1Ref, xInlierVec2Ref;
        Vector2dVec pInlierVec2Ref;
        vector<double> sizeVec2Ref;
        
        // select inliers
        for (int i = 0; i < reprojectedVec.size(); i++)
        {
            if (errVec[i] < 3.6*sigmaSq)
            {
                xInlierVec1Ref.push_back(xInlierVec1[i]);
                xInlierVec2Ref.push_back(xInlierVec2[i]);
                pInlierVec2Ref.push_back(pInlierVec2[i]);
                sizeVec2Ref.push_back(sizeVec2[i]);
            }
        }
        
        computeTransfSparse(xInlierVec1Ref, xInlierVec2Ref, pInlierVec2Ref, sizeVec2Ref,
                xiOdom.inverseCompose(xiOdomNew), xiIncr, true);
        
        xiLocal = xiLocal.compose(xiIncr);
        
        if (abs(xiOdom.inverseCompose(xiOdomNew).trans().norm() - xiIncr.trans().norm()) / xiIncr.trans().norm() > 0.05 or 0) //to display the RANSAC result
        {
//            Mat img_matches;
//            drawMatches( imageOld, keypointVec1, imageNew, keypointVec2,
//               goodMatchVec, img_matches, Scalar::all(-1), Scalar::all(-1),
//               vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
//            imshow("matches", img_matches);

            Mat8uc3 lineMat;
             cvtColor(imageNew, lineMat, CV_GRAY2BGR);
            Transf dxi = xiBaseCam.inverseCompose(xiIncr).compose(xiBaseCam);
            Vector2dVec reprojectedVec = reprojectPoints(reconstVec1, reconstVec2, dxi);
            
            for (auto & x : pInlierVec2Ref)
            {
                cross(lineMat, x[0], x[1], 3, Scalar(0, 0, 255), 3);
            }
            
            for (int i = 0; i < inlierMask.size(); i++)
            {   
                if (inlierMask[i]) 
                {
                    
                    cross(lineMat, reprojectedVec[i][0], reprojectedVec[i][1], 3, Scalar(250, 0, 0), 1);
                    line(lineMat, Point(pointVec1[i][0], pointVec1[i][1]),
                                        Point(pointVec2[i][0], pointVec2[i][1]), Scalar(0, 255, 0));
                    
//                    cout << setw(5) << i << setw(15) << pointVec2[i].transpose() 
//                            << setw(15) << reprojectedVec[i].transpose() << setw(15) << sizeVec[i] << endl;
//                    if ((pointVec1[i] - pointVec2[i]).norm() > 25)
//                    {
//                        imshow("lines", lineMat);
//                        waitKey();
//                    }
                }
                
            }
            imshow("lines", lineMat);
            waitKey(1500);
        }
        
    }
    
    //refresh state
    keypointVec1.swap(keypointVec2);
    swap(desc1, desc2);
    imageNew.copyTo(imageOld);
    xiOdom = xiOdomNew;
}
    
      
double SparseOdometry::computeTransfSparse(const Vector3dVec & xVec1, const Vector3dVec & xVec2, 
        const Vector2dVec & pVec2, const vector<double> & sizeVec, const Transf xiOdom, Transf & xiOut, bool report)
{
    Problem problem;
    //TODO avoid reconstruction in first place
    
    array<double, 6> xiArr = xiOdom.toArray();
    
    OdometryPrior * odometryCost = new OdometryPrior(0.03, 0.5, 0.03, 0.05, xiOdom);
    problem.AddResidualBlock(odometryCost, NULL, xiArr.data());
    
    //TODO use inlier mask instead
    SparseReprojectCost * projectionCost = new SparseReprojectCost(camera,
                                                    xVec1, xVec2, pVec2, sizeVec, xiBaseCam);
    problem.AddResidualBlock(projectionCost, NULL, xiArr.data());
    
    Solver::Options options;
//        options.check_gradients = true;
    Solver::Summary summary;
//        options.function_tolerance = 1e-8;
//        options.gradient_tolerance = 1e-8;
//        options.parameter_tolerance = 1e-8;
    options.max_num_iterations = 25;
//        options.numeric_derivative_relative_step_size = 1e-8;
    Solve(options, &problem, &summary);
//        cout << summary.BriefReport() << endl;
//    if (report) cout << summary.FullReport() << endl;
    xiOut = Transf(xiArr.data());
    return summary.final_cost; 
}


Vector2dVec SparseOdometry::reprojectPoints(const Vector3dVec & cloud1,
    const Vector3dVec & cloud2, const Transf dxi)
{    
    vector<double> lambdaVec(cloud1.size());
    Triangulator(dxi).computeRegular(cloud1, cloud2, lambdaVec.data());
    Vector3dVec cloudScaled;
    cloudScaled.reserve(cloud1.size());
//    cout << "LAMBDAS" << endl;
    for (int i = 0; i < cloud1.size(); i++)
    {
//        cout << setw(5) << i << setw(15) << lambdaVec[i] << endl;
        cloudScaled.emplace_back(cloud1[i] * lambdaVec[i]);
    }
    
    Vector3dVec cloudTransformed;
    
    dxi.inverseTransform(cloudScaled, cloudTransformed);
//    cout << "TRANSFORMED" << endl;
//    for (int i = 0; i < cloud1.size(); i++)
//    {
//        cout    << setw(5) << i << setw(15) << cloudScaled[i].transpose() 
//                << setw(15) << cloudTransformed[i].transpose() << endl;
//    }
    
    Vector2dVec ptCheckVec;
    camera->projectPointCloud(cloudTransformed, ptCheckVec);
    
//    Vector2dVec ptOrigVec; //FIXME debug
//    camera->projectPointCloud(cloud2, ptOrigVec);
//    for (int i = 0; i < cloud1.size(); i++)
//    {
//        cout    << setw(5) << i << setw(15) << ptOrigVec[i].transpose() 
//                << setw(15) << ptCheckVec[i].transpose() << endl;
//    }
    
    return ptCheckVec;
}

   
void SparseOdometry::ransacNPoints(const Vector3dVec & cloud1,
    const Vector3dVec & cloud2, const Vector2dVec & ptVec2, const vector<double> & sizeVec,
    const Transf xiOdom, vector<bool> & inlierMask)
{
    assert(cloud1.size() == cloud2.size());
    //define constants
    const int maxIteration = 200;
    const double thresh = 1.;
    int inlierCount = numRansacPoints;
    
    vector<int> indexVec;
    for (int idx = 0; idx < cloud1.size(); idx++)
    {
        indexVec.push_back(idx);
    }
    
    //triangulate all the vectors
    Transf dxi = xiBaseCam.inverseCompose(xiOdom).compose(xiBaseCam);
    
    
    for (int ransacIteration = 0; ransacIteration < maxIteration; ransacIteration++)
    {
        //TODO optimiza the step
        shuffle(indexVec.begin(), indexVec.end(), _g);
        Vector3dVec xSampleVec1, xSampleVec2;
        Vector2dVec pSampleVec2;
        vector<double> sizeVec2;
        for (auto iter = indexVec.begin(); iter != indexVec.begin() + numRansacPoints; ++iter)
        {
            xSampleVec1.push_back(cloud1[*iter]);
            xSampleVec2.push_back(cloud2[*iter]);
            pSampleVec2.push_back(ptVec2[*iter]);
            sizeVec2.push_back(sizeVec[*iter]);
//            sizeVec2.push_back(1);
        }
        
        // fit the model
        Transf xiOut;
        //TODO // chi2 test, 16 variables, 2% confidence
        computeTransfSparse(xSampleVec1, xSampleVec2, pSampleVec2, sizeVec2, xiOdom, xiOut);
        
        vector<double> lambdaVec(cloud1.size());
        xiOut = xiBaseCam.inverseCompose(xiOut).compose(xiBaseCam);
        Triangulator(xiOut).computeRegular(cloud1, cloud2, lambdaVec.data());
        Vector3dVec cloudScaled;
        cloudScaled.reserve(cloud1.size());
        for (int i = 0; i < cloud1.size(); i++)
        {
            cloudScaled.emplace_back(cloud1[i] * lambdaVec[i]);
        }
        xiOut.inverseTransform(cloudScaled, cloudScaled);
        Vector2dVec ptCheckVec;
        camera->projectPointCloud(cloudScaled, ptCheckVec);
        // count inliers
        vector<bool> inlierMaskEstim;
        int inlierCountEstim = 0;
        
        for (int idx = 0; idx < ptCheckVec.size(); idx++)
        {
            double res = (ptVec2[idx] - ptCheckVec[idx]).norm();
            if (res < thresh)
            {
                inlierCountEstim++;
                inlierMaskEstim.push_back(true);
            }
            else
            {
                inlierMaskEstim.push_back(false);
            }
        }
        
        // refresh the best hypothesis
        if (inlierCountEstim > inlierCount)
        {
            
            
            //FIXME DEBUG
            cout << "RANSAC IMPROVE  " << ransacIteration << setw(10) << inlierCountEstim << endl;
//            
//            for (int i = 0; i < ptCheckVec.size(); i++)
//            {
//                if (inlierMaskEstim[i])
//                {
//                    cout << setw(5) << i << setw(15) << ptVec2[i].transpose() 
//                            << setw(15) << ptCheckVec[i].transpose() << endl; 
//                }
//            }
            
            inlierCount = inlierCountEstim;
            inlierMask.swap(inlierMaskEstim);
//            if (inlierCount > 30) break;
        }
    }
//    cout << "inlierCount : " << inlierCount << endl;
}


            
