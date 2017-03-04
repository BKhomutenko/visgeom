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
    
#include "calibration/corner_detector.h"

#include "eigen.h"
#include "ceres.h"
#include "ocv.h"
#include "io.h"

#include "utils/curve_rasterizer.h"

double findMinDistance(const Vector2dVec & cornerVec, const int rows, const int cols)
{
    assert(cornerVec.size() == rows * cols);
    double res = 1e10;
    //horizontal distance
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols - 1; j++)
        {
            cout << i << ' ' << j << endl;
            int idx = i * rows + i;
            res = min(res, (cornerVec[idx] - cornerVec[idx + 1]).squaredNorm());
        }
    }
    
    for (int i = 0; i < rows - 1; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int idx = i * rows + i;
            res = min(res, (cornerVec[idx] - cornerVec[idx + cols]).squaredNorm());
        }
    }
    return sqrt(res);
}

void CornerDetector::improveCorners(Vector2dVec & pointVec) const
{
    for (auto & x : pointVec)
    {
        array<double, 4> dataArr;
        initPoin(x, dataArr.data());
        
        ceres::GradientProblem problem(new SubpixelCorner(_gradx, _grady));
        ceres::GradientProblemSolver::Options options;
//        options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
//        options.line_search_type = ceres::ARMIJO;
//        options.line_search_sufficient_function_decrease = 1e-6;
//        options.max_num_line_search_step_size_iterations = 5;
        options.logging_type = ceres::SILENT;
//        options.minimizer_progress_to_stdout = true;
        ceres::GradientProblemSolver::Summary summary;
        ceres::Solve(options, problem, dataArr.data(), &summary);
        
//        cout << summary.FullReport() << endl;
        
//        cout << dataArr[0] << "   "
//            << dataArr[1] << "   "
//            << dataArr[2] << "   "
//            << dataArr[3] << endl;
            
        x[0] = dataArr[0];
        x[1] = dataArr[1];
    }
}

CornerDetector::CornerDetector(const Mat8u & img, const int initRadius) : //: _img(img)
    _gradx(img.rows, img.cols, CV_32F),
    _grady(img.rows, img.cols, CV_32F),
    _img(img.rows, img.cols, CV_32F),
    INIT_RADIUS(initRadius)
{
    Sobel(img, _gradx, CV_32F, 1, 0, 3, 1e-3);
    Sobel(img, _grady, CV_32F, 0, 1, 3, 1e-3);
    img.copyTo(_img);
//    Mat32f Ixx(img.rows, img.cols, CV_32F), Ixy(img.rows, img.cols, CV_32F);
//    Mat32f Iyx(img.rows, img.cols, CV_32F), Iyy(img.rows, img.cols, CV_32F);
//    
//    img.convertTo(_img, CV_32F);
//    Grid2D imgGrid(_img.cols, _img.rows, (float*)_img.data);
//    ceres::BiCubicInterpolator<Grid2D> imgInter(imgGrid);
//    for (int r = 0; r < img.rows; r++)
//    {
//        for (int c = 0; c < img.cols; c++)
//        {
//            double foo, fx, fy;
//            imgInter.Evaluate(r, c, &foo, &fy, &fx);
//            _gradx(r, c) = fx / 10;
//            _grady(r, c) = fy / 10;
//        }
//    }
//    
//    const double delta = 0.5;
//    
//    //Ixx Ixy
//    Grid2D gxGrid(_img.cols, _img.rows, (float*)_gradx.data);
//    ceres::BiCubicInterpolator<Grid2D> gxInter(gxGrid);
//    for (int r = 0; r < img.rows; r++)
//    {
//        for (int c = 0; c < img.cols; c++)
//        {
//            double foo, fx, fy;
//            gxInter.Evaluate(r + delta, c + delta, &foo, &fy, &fx); 
//            Ixy(r, c) = fy+0.5;
//            Ixx(r, c) = fx+0.5;
//        }
//    }
//    
//    //Iyx Iyy
//    Grid2D gyGrid(_img.cols, _img.rows, (float*)_grady.data);
//    ceres::BiCubicInterpolator<Grid2D> gyInter(gyGrid);
//    for (int r = 0; r < img.rows; r++)
//    {
//        for (int c = 0; c < img.cols; c++)
//        {
//            double foo, fx, fy;
//            gyInter.Evaluate(r + delta, c + delta, &foo, &fy, &fx); 
//            Iyy(r, c) = fy+0.5;
//            Iyx(r, c) = fx+0.5;
//        }
//    }
//    
//    imshow("img", img);
//    imshow("gradx", _gradx);
//    imshow("grady", _grady);
//    
//    imshow("Ixy", Ixy);
//    imshow("Iyx", Iyx);
//    imshow("err", Iyx - Ixy + 0.5);
//    
//    waitKey();
}

void CornerDetector::initPoin(const Vector2d & pt, double * data) const
{
    // for init_radius look for transitions
    
    Polynomial2 circle;
    
    circle.kuu = 1;
    circle.kvv = 1;
    circle.kuv = 0;
    circle.ku = -2*pt[0];
    circle.kv = -2*pt[1];
    circle.k1 = pt[0] * pt[0] + pt[1] * pt[1] - INIT_RADIUS * INIT_RADIUS;
    
    Vector2i pti = round(pt);
    CurveRasterizer<int, Polynomial2> raster(pti[0] + INIT_RADIUS, pti[1],
                                             pti[0], pti[1] + INIT_RADIUS, circle);
    //////////
//    Mat32f img(_img.size());
//    
//    _img.copyTo(img);
//    cout << img.size() << endl;

    //go along the circle and save the transition strengths
    vector<double> transitionVec;
    Vector2iVec circleVec;
    for (int i = 0; i < ceil(2 * INIT_RADIUS * M_PI); i++, raster.step())
    {
        Vector2d r(raster.u - pti[0], raster.v - pti[1]);
        Vector2d g(_gradx(raster.v, raster.u), _grady(raster.v, raster.u));
        transitionVec.push_back(g[0] * r[1] - g[1] * r[0]);
        circleVec.emplace_back(raster.u, raster.v);
        
        /////////
//        img(raster.v, raster.u) = 0;

    }
    
    ////find the max elements
    auto maxIter1 = max_element(transitionVec.begin(), transitionVec.end());
    
    //find the max element on the other side of the circle
    int maxIdx1 = distance(transitionVec.begin(), maxIter1);
    int searchLimit1 = maxIdx1 + transitionVec.size() / 4;
    int searchLimit2 = searchLimit1 + transitionVec.size() / 2;
    
    //TODO optimize in future
    int maxIdx2;
    if (searchLimit1 < transitionVec.size() and searchLimit2 >= transitionVec.size())
    {
        searchLimit2 = searchLimit2 % transitionVec.size();
        auto maxIter21 = max_element(transitionVec.begin() + searchLimit1,
                                    transitionVec.end());
        auto maxIter22 = max_element(transitionVec.begin(),
                                    transitionVec.begin() + searchLimit2); 
        if (*maxIter21 > *maxIter22) maxIdx2 = distance(transitionVec.begin(), maxIter21);
        else maxIdx2 = distance(transitionVec.begin(), maxIter22);
    } 
    else
    {
        searchLimit1 = searchLimit1 % transitionVec.size();
        searchLimit2 = searchLimit2 % transitionVec.size();
        auto maxIter2 = max_element(transitionVec.begin() + searchLimit1,
                                    transitionVec.begin() + searchLimit2);
        maxIdx2 = distance(transitionVec.begin(), maxIter2);
    }
    
    ////find the min elements
    
    auto minIter1 = min_element(transitionVec.begin(), transitionVec.end());
    
    //find the min element on the other side of the circle
    int minIdx1 = distance(transitionVec.begin(), minIter1);
    searchLimit1 = minIdx1 + transitionVec.size() / 4;
    searchLimit2 = searchLimit1 + transitionVec.size() / 2;
    
    //TODO optimize in future
    int minIdx2;
    if (searchLimit1 < transitionVec.size() and searchLimit2 >= transitionVec.size())
    {
        searchLimit2 = searchLimit2 % transitionVec.size();
        auto minIter21 = min_element(transitionVec.begin() + searchLimit1,
                                    transitionVec.end());
        auto minIter22 = min_element(transitionVec.begin(),
                                    transitionVec.begin() + searchLimit2); 
        if (*minIter21 < *minIter22) minIdx2 = distance(transitionVec.begin(), minIter21);
        else minIdx2 = distance(transitionVec.begin(), minIter22);
    } 
    else
    {
        searchLimit1 = searchLimit1 % transitionVec.size();
        searchLimit2 = searchLimit2 % transitionVec.size();
        auto minIter2 = min_element(transitionVec.begin() + searchLimit1,
                                    transitionVec.begin() + searchLimit2);
        minIdx2 = distance(transitionVec.begin(), minIter2);
    }
    
    
    
    //compute the angles and the intersection of the lines
    Vector2i A = circleVec[maxIdx1];
    Vector2i B = circleVec[minIdx1];
    Vector2i C = circleVec[maxIdx2];
    Vector2i D = circleVec[minIdx2];
    
    ///////////////
//    img(A[1], A[0]) = 255;
//    img(B[1], B[0]) = 255;
//    img(C[1], C[0]) = 255;
//    img(D[1], D[0]) = 255;
//    
//    imshow("transitions", img / 255);
//    waitKey();

    //to find the intersection E we solve:
    // {AC x AE = 0
    // {BD x BE = 0 
    Matrix2d M;
    Vector2d b;
    M(0, 0) = A[1] - C[1];
    M(0, 1) = C[0] - A[0];
    M(1, 0) = B[1] - D[1];
    M(1, 1) = D[0] - B[0];
    b[0] = A[1] * M(0, 1) + A[0] * M(0, 0);
    b[1] = B[1] * M(1, 1) + B[0] * M(1, 0);
    Vector2d E = M.inverse() * b;
    data[0] = E[0];
    data[1] = E[1];
    data[2] = atan2(A[1] - C[1], A[0] - C[0]);
    data[3] = atan2(B[1] - D[1], B[0] - D[0]);
}
    
//    Mat32f _gradx, _grady;
//    Mat8u _img;
    
