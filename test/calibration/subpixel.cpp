/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your opton) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#include "io.h"
#include "ocv.h"
#include "eigen.h"

//#include "calibration/corner_detector.h"

#include "utils/curve_rasterizer.h"

template<typename T>
void setZero(T begin, T end, T it)
{
    *it = 0;
    if (it == begin) *(end - 1) = 0;
    else *(it - 1) = 0;
    
    if (it + 1 == end) *begin = 0;
    else *(it + 1) = 0;
}

bool comp(const pair<double, Vector2i> & a, const pair<double, Vector2i> & b)
{
    return a.first < b.first;
}

int main(int argc, char** argv) {

    const int Nx = 8;
    const int Ny = 5;
    
    vector<cv::Point2f> centers;
    Size patternSize(Nx, Ny);
    Vector2dVec cornerVec;
    Mat8u frame = imread(argv[1], 0);
    
    imshow("frame", frame);
    
//    resize(frame, frame, Size(0, 0), 2, 2);
    
    ///TEST A RESPONSE FUNCTION
    
    Mat8u src;
    const double SIGMA = 2.5;
    const int FILTER_SIZE = 1 + 2 * round(SIGMA);
    GaussianBlur(frame, src, Size(FILTER_SIZE, FILTER_SIZE), SIGMA, SIGMA);
    
    Mat32f respXY(src.size()), respUV(src.size()), resp(src.size());
    respXY.setTo(0);
    respUV.setTo(0);
    resp.setTo(0);
    double acc = 0;
    int count = 0;
    const int SIZE = 3;
    const int HSIZE = 2;
    for (int v = SIZE; v < src.rows - SIZE; v++)
    {
        for (int u = SIZE; u < src.cols - SIZE; u++)
        {
            
            double gx = (src(v, u + HSIZE) - src(v, u - HSIZE) ) / (2 * HSIZE);
            double gy = (src(v + HSIZE, u) - src(v - HSIZE, u) ) / (2 * HSIZE);
            
            
            //We are looking for saddle points, that is negative hessian determinant
            double iuu = src(v, u - SIZE) + src(v, u + SIZE) - 2* src(v, u);
            iuu /= SIZE*SIZE;
            double ivv = src(v - SIZE, u) + src(v + SIZE, u) - 2* src(v, u);
            ivv /= SIZE*SIZE;
            double iuv = src(v - HSIZE, u - HSIZE) + src(v + HSIZE, u + HSIZE) 
                            - src(v  + HSIZE, u - HSIZE) - src(v - HSIZE, u + HSIZE);
            iuv /= 4*HSIZE*HSIZE;
            respXY(v, u) = -iuu * ivv;
            respUV(v, u) = iuv * iuv ;
            double respVal = -iuu * ivv + iuv * iuv - 0.01*pow(gx * gx + gy * gy, 2);
            if (respVal > 1) 
            {
                resp(v, u) = respVal;
                acc += respVal;
                count++;
            }
        }
    }
    const double avgVal = acc / count;
    imshow("respXY", -respXY / 25);
    imshow("respUV", respUV / 25);
    imshow("resp", resp / 250);
    
    waitKey();
    
    
    
    
    ///TEST HARRIS CORNERS
    const int IMAGE_MARGIN = 2;
    Mat32f gradx, grady, respxx, respxy, respyy;
    
    Sobel(frame, gradx, CV_32F, 1, 0, 3, 1./255);
    Sobel(frame, grady, CV_32F, 0, 1, 3, 1./255);
    respxx = gradx.mul(gradx);
    respyy = grady.mul(grady);
    respxy = gradx.mul(grady);
    
    const double SIGMA2 = 0.8;
    GaussianBlur(respxx, respxx, Size(5, 5), SIGMA2, SIGMA2);
    GaussianBlur(respxy, respxy, Size(5, 5), SIGMA2, SIGMA2);
    GaussianBlur(respyy, respyy, Size(5, 5), SIGMA2, SIGMA2);
    
//    Mat32f resp(frame.size());
//    resp.setTo(0);
//    const double KAPPA = 0.05;
//    double acc = 0;
//    int count = 0;
//    for (int v = IMAGE_MARGIN; v < frame.rows - IMAGE_MARGIN; v++)
//    {
//        for (int u = IMAGE_MARGIN; u < frame.cols - IMAGE_MARGIN; u++)
//        {
//            const double ixx = respxx(v, u);
//            const double ixy = respxy(v, u);
//            const double iyy = respyy(v, u);
//            const double respVal = ixx * iyy - ixy * ixy - KAPPA * pow(ixx + iyy, 2);
//            if (respVal > 0)
//            {
//                resp(v, u) = respVal;
//                acc += respVal*respVal;
//                count++;
//            }
//        }
//        
//    }
//    const double avgVal = sqrt(acc / count);
//    cout << avgVal << endl;
//    imshow("resp2", resp);
//    waitKey();
    
    //find maxima
    
    
    Mat32f detected(frame.size());
    detected.setTo(0);
    
    const int WINDOW_SIZE = 4;
    vector<pair<double, Vector2i>> hypVec;
    for (int v = WINDOW_SIZE; v < frame.rows - WINDOW_SIZE; v++)
    {
        for (int u = WINDOW_SIZE; u < frame.cols - WINDOW_SIZE; u++)
        {
            
            bool isMax = true;
            const double val = resp(v, u);
            if (val < avgVal) continue;
            for (int j = -WINDOW_SIZE; j <= WINDOW_SIZE; j++)
            {
                for (int i = -WINDOW_SIZE; i <= WINDOW_SIZE; i++)
                {
                    if (i == 0 and j == 0) continue;
                    const double valNeigh = resp(v + j, u + i);
                    if ( val <= valNeigh)
                    {
                        //quite unlikely scenario that 
                        //the two neighbor pixels are local maxima simultaneously
                        //in this case only one is to be selected
                        if (val == valNeigh and (i > 0 or i == 0 and j > 0)) continue;
                        isMax = false;
                        break;
                    }
                    
                }
                if (not isMax) break;
            }
            if (isMax)
            { 
                hypVec.emplace_back(val, Vector2i(u, v));
//                detected(v, u) = 1;
            }
        }
    }
    cout << hypVec.size() << endl;
    
    make_heap(hypVec.begin(), hypVec.end(), comp);
    /// Check coner-ness
   
    
    
    const int INIT_RADIUS = 5;
    
    
                                              
    
    int count2 = 0;
    
    for (int i = 0; i < 60; i++)
    {
        if (hypVec.empty()) break;
        pop_heap(hypVec.begin(), hypVec.end(), comp);
        auto pt = hypVec.back().second;
        hypVec.pop_back();
        // compute transitions
        Polynomial2 circle;
        circle.kuu = 1;
        circle.kvv = 1;
        circle.kuv = 0;
        circle.ku = -2*pt[0];
        circle.kv = -2*pt[1];
        circle.k1 = pt[0] * pt[0] + pt[1] * pt[1] - INIT_RADIUS * INIT_RADIUS;
        
        Vector2i pt0(pt[0] + INIT_RADIUS, pt[1]);
        CurveRasterizer<int, Polynomial2> raster(pt[0] + INIT_RADIUS, pt[1],
                                                 pt[0], pt[1] + INIT_RADIUS, circle);
        
        vector<double> sampleVec;
        Vector2iVec circleVec;
        for (int i = 0;
             not (i > 5 and abs(raster.u - pt0[0]) <= 1 and abs(raster.v - pt0[1]) <= 1);
             i++, raster.step())
        {
            sampleVec.push_back(frame(raster.v, raster.u));
            circleVec.emplace_back(raster.u, raster.v);
        }
        
        double totalVal = 0;
        double diffVal = 0;
        
//        int deltaStep = sampleVec.size() / 2;
//        for (int i = 0; i < sampleVec.size(); i++)
//        {
//            int j = (i  + deltaStep) % sampleVec.size();
//            totalVal += sampleVec[i];
//            diffVal += abs(sampleVec[i] - sampleVec[j]);
//        }
//        if (diffVal / totalVal > 0.2) continue;
//        detected(pt[1], pt[0]) = 1;
//        count2++;
                                 
//        vector<double> transitionVec;
//        
//        transitionVec.push_back(sampleVec[1] - sampleVec.back());
//        
//        for (int i = 1; i < sampleVec.size() - 1; i++)
//        {
//            transitionVec.push_back(sampleVec[i + 1] - sampleVec[i - 1]);
//        }
//        int idxLast = sampleVec.size() - 2;
//        transitionVec.push_back(sampleVec.front() - sampleVec[idxLast]);
//                                
//        // find the strongest
//        
//        auto itMax = max_element(transitionVec.begin(), transitionVec.end());
//        double transMax = *itMax;
//        
//        // remove the maximum and its neighbors
//        setZero(transitionVec.begin(), transitionVec.end(), itMax);
//        // check the other three : must be at least 0.75
//        
//        auto itMax2 = max_element(transitionVec.begin(), transitionVec.end());
//        
//        if (*itMax2 < transMax * 0.75) continue;
//        setZero(transitionVec.begin(), transitionVec.end(), itMax2);
        
//        auto itMax3 = max_element(transitionVec.begin(), transitionVec.end());
//        
//        if (*itMax3 > transMax * 0.3) continue;
//        setZero(transitionVec.begin(), transitionVec.end(), itMax3);
        
//        auto itMin = min_element(transitionVec.begin(), transitionVec.end());
//        
//        if (*itMin > transMax * -0.75) continue;
//        setZero(transitionVec.begin(), transitionVec.end(), itMin);
//        
//        auto itMin2 = min_element(transitionVec.begin(), transitionVec.end());
//        
//        if (*itMin2 > transMax * -0.75) continue;
//        setZero(transitionVec.begin(), transitionVec.end(), itMin2);
        
//        auto itMin3 = min_element(transitionVec.begin(), transitionVec.end());
//        
//        if (*itMin3 < transMax * -0.3) continue;
//        setZero(transitionVec.begin(), transitionVec.end(), itMin3);
        
        detected(pt[1], pt[0]) = 1;
        count2++;
    }
    cout << count2 << endl;
    imshow("detected", detected);
    waitKey();
    return 0;
    
    
    
    /////////////////////////////////////////////////////////////////////
//    Mat8u corners1, corners2;
//    frame.copyTo(corners1);
//    frame.copyTo(corners2);
//    bool patternIsFound = findChessboardCorners(frame, patternSize, centers, CV_CALIB_CB_ADAptVE_THRESH);
//    if (not patternIsFound)
//    {
//        cout << argv[1] << " : ERROR, pattern not found" << endl;
//        return 0;
//    }

//    
//    
//        
//    cornerVec.resize(Nx * Ny);
//    for (int i = 0; i < Nx * Ny; i++)
//    {
//        cornerVec[i] = Vector2d(centers[i].x, centers[i].y);
//    }
//    
//    double minDist = findMinDistance(cornerVec, Ny, Nx);
//    
//    CornerDetector detector(frame, min(minDist / 3., 10.));
//    
//    cout << "before : " << endl;
//    for (auto & pt : cornerVec)
//    {
//        cout << pt.transpose() << endl;
//    }
//    
//    
//    
//    detector.improveCorners(cornerVec);
//    
//    cout << "after : " << endl;
//    for (auto & pt : cornerVec)
//    {
//        cout << pt.transpose() << endl;
//    }
//    
//    drawChessboardCorners(corners1, patternSize, Mat(centers), patternIsFound);
//    imshow("corners1", corners1);
//    
//    vector<cv::Point2f> centers2;
//    for (auto & x : cornerVec)
//    {
//        centers2.emplace_back(x[0], x[1]);
//    }
//    
//    drawChessboardCorners(corners2, patternSize, Mat(centers2), patternIsFound);
//    imshow("corners2", corners2);
//    
//    waitKey();
//    return 0;
}
