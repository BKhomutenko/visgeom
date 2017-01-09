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
#include "eigen.h"
#include "ocv.h"
#include "ceres.h" 

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "reconstruction/triangulator.h"
#include "localization/local_cost_functions.h"
    
void SparseOdometry::feedData(const Mat8u & imageNew, const Transf xiOdomNew)
{
    Transf dxi = xiBaseCam.inverseCompose(xiOdom.inverseCompose(xiOdomNew)).compose(xiBaseCam);
    if (dxi.trans().norm() < MIN_STEREO_BASE and not keypointVec1.empty()) return;
    
    keypointVec2.clear();
    detector.detect(imageNew, keypointVec2);
    detector.compute(imageNew, keypointVec2, desc2);
    cout << "detected : " << keypointVec2.size() << endl; 
    if (not keypointVec1.empty())
    {
        //correspondence
        BFMatcher matcher(cv::NORM_HAMMING, true);
        vector<DMatch> matchVec;
        matcher.match(desc1, desc2, matchVec);
        
        //reconstruct
        vector<DMatch> goodMatchVec;
        Vector2dVec pointVec1, pointVec2;
        for (auto & match : matchVec)
        {
            if (match.distance > distThresh) continue;
            
            goodMatchVec.push_back(match);
            auto pt1 = keypointVec1[match.queryIdx].pt;
            auto pt2 = keypointVec2[match.trainIdx].pt;
            
//                if (norm(pt1 - pt2) < 3) continue;
            
            pointVec1.emplace_back(pt1.x, pt1.y);
            pointVec2.emplace_back(pt2.x, pt2.y);
        }
        cout << endl;
        Vector3dVec reconstVec1, reconstVec2;
        camera->reconstructPointCloud(pointVec1, reconstVec1);
        camera->reconstructPointCloud(pointVec2, reconstVec2);
        
        //RANSAC
        vector<bool> inlierMask;
        ransacFivePoint(reconstVec1, reconstVec2, pointVec2,
                xiOdom.inverseCompose(xiOdomNew), inlierMask);
        
        if (1) //to display the RANSAC result
        {
            Mat img_matches;
            drawMatches( imageNew, keypointVec1, imageNew, keypointVec2,
               goodMatchVec, img_matches, Scalar::all(-1), Scalar::all(-1),
               vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
            
            Mat8u lineMat;
            imageNew.copyTo(lineMat);
            lineMat *= 0.8;
            
            imshow("matches", img_matches);
            for (int i = 0; i < inlierMask.size(); i++)
            {   
                if (inlierMask[i]) 
                {
                    line(lineMat, Point(pointVec1[i][0], pointVec1[i][1]),
                                        Point(pointVec2[i][0], pointVec2[i][1]), 255);
                }
            }
            imshow("lines", lineMat);
            
        }
        
        
        //refinement
        Vector3dVec xInlierVec1, xInlierVec2;
        Vector2dVec pInlierVec2;
        for (int i = 0; i < inlierMask.size(); i++)
        {
            if (not inlierMask[i]) continue;
            
            xInlierVec1.push_back(reconstVec1[i]);
            xInlierVec2.push_back(reconstVec2[i]);
            pInlierVec2.push_back(pointVec2[i]);
        }
        Transf xiOut;
        computeTransfSparse(xInlierVec1, xInlierVec2, pInlierVec2, 
                xiOdom.inverseCompose(xiOdomNew), xiOut, true);
        cout << xiOut << endl;
        cout << xiOdom.compose(xiOut) << endl;
        waitKey();
    }
    
    //refresh state
    keypointVec1.swap(keypointVec2);
    swap(desc1, desc2);
    xiOdom = xiOdomNew;
}
    
      
double SparseOdometry::computeTransfSparse(const Vector3dVec & xVec1, const Vector3dVec & xVec2, 
        const Vector2dVec & pVec2, const Transf xiOdom, Transf & xiOut, bool report)
{
    Problem problem;
    //TODO avoid reconstruction in first place
    
    array<double, 6> xiArr = xiOdom.toArray();
    
    OdometryPrior * odometryCost = new OdometryPrior(0.03, 0.03, 0.01, 0.03, xiOdom);
    problem.AddResidualBlock(odometryCost, NULL, xiArr.data());
    
    SparseReprojectCost * projectionCost = new SparseReprojectCost(camera,
                                                    xVec1, xVec2, pVec2, xiBaseCam);
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
    if (report) cout << summary.FullReport() << endl;
    xiOut = Transf(xiArr.data());
    return summary.final_cost; 
}
    
void SparseOdometry::ransacFivePoint(const Vector3dVec & cloud1,
    const Vector3dVec & cloud2, const Vector2dVec & ptVec2,
    const Transf xiOdom, vector<bool> & inlierMask)
{
    assert(cloud1.size() == cloud2.size());
    //define constants
    int maxIteration = 50;
    double thresh = 3;
    const int numRansacPoints = 2;
    int inlierCount = numRansacPoints;
    
    vector<int> indexVec;
    for (int idx = 0; idx < cloud1.size(); idx++)
    {
        indexVec.push_back(idx);
    }
    
    //triangulate all the vectors
    Transf dxi = xiBaseCam.inverseCompose(xiOdom).compose(xiBaseCam);
    
    mt19937 g(0);
    for (int ransacIteration = 0; ransacIteration < maxIteration; ransacIteration++)
    {
        //TODO optimiza the step
        shuffle(indexVec.begin(), indexVec.end(), g);
        Vector3dVec xSampleVec1, xSampleVec2;
        Vector2dVec pSampleVec2;
        for (auto iter = indexVec.begin(); iter != indexVec.begin() + numRansacPoints; ++iter)
        {
            xSampleVec1.push_back(cloud1[*iter]);
            xSampleVec2.push_back(cloud2[*iter]);
            pSampleVec2.push_back(ptVec2[*iter]);
        }
        
        // fit the model
        Transf xiOut;
        //TODO // chi2 test, 16 variables, 2% confidence
        computeTransfSparse(xSampleVec1, xSampleVec2, pSampleVec2, xiOdom, xiOut);
        
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
            double res = (ptVec2[idx] - ptCheckVec[idx]).squaredNorm();
            if (abs(res) < thresh)
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
            inlierCount = inlierCountEstim;
            inlierMask.swap(inlierMaskEstim);
//            if (inlierCount > 30) break;
        }
    }
    cout << "inlierCount : " << inlierCount << endl;
}


            
