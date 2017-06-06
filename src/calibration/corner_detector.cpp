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

//SubpixelCorner function definition

SubpixelCorner::SubpixelCorner(const Mat32f & gradu, const Mat32f & gradv, const int steps, const double length) :    
            _graduGrid(gradu.cols, gradu.rows, (float*)gradu.data),
            _gradvGrid(gradv.cols, gradv.rows, (float*)gradv.data),
            _stepLength(length / steps)
{
    stepVec.reserve(2*steps);
    for (int i = 1; i <= steps; i++)
    {
        stepVec.push_back(-i * _stepLength);
        stepVec.push_back(i * _stepLength);
    }
}
            
bool SubpixelCorner::Evaluate(const double* parameters,
                    double* cost,
                    double* gradient) const 
{
    const double & u = parameters[0];
    const double & v = parameters[1];
    
    //TODO possibly merge
    *cost = 0;
    if (gradient != NULL) fill(gradient, gradient + 5, 0.);
    
    ceres::BiCubicInterpolator<Grid2D<float>> graduInter(_graduGrid);
    ceres::BiCubicInterpolator<Grid2D<float>> gradvInter(_gradvGrid);
    for (int direction = 0; direction < 2; direction++)
    {
        int thIdx = 2 + direction;
        const double s = sin(parameters[thIdx]);
        const double c = cos(parameters[thIdx]);
        const double & h = parameters[4];
        double flowDir = direction ? 1 : -1;
        
        for (double length : stepVec)
        {
            double eta = (length > 0 ? 1 : -1) * flowDir;
            double ui = u + c * length - s * h * eta;
            double vi = v + s * length + c * h * eta;
            
            double fu, fuu, fuv;
            graduInter.Evaluate(vi , ui,
                        &fu, &fuv, &fuu);
            double fv, fvu, fvv;
            gradvInter.Evaluate(vi , ui,
                        &fv, &fvv, &fvu);
            *cost += eta*(fv * c - fu * s);
            if (gradient != NULL)
            {
                double dudth = -s * length - c * h * eta;
                double dvdth = c * length - s * h * eta;
                gradient[0] += eta * (fvu * c - fuu * s);
                gradient[1] += eta * (fvv * c - fuv * s);
                gradient[thIdx] += eta * ( (fvv * dvdth + fvu * dudth) * c 
                                - (fuv * dvdth  + fuu * dudth) * s
                                - fu * c - fv * s );
                gradient[4] += fvv * c * c + fuu * s * s - s * c * (fvu + fuv); 
            }
        }
    }
    return true;
}




double findMinDistance(const Vector2dVec & cornerVec, const int rows, const int cols)
{
    assert(cornerVec.size() == rows * cols);
    double res = 1e10;
    //horizontal distance
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols - 1; j++)
        {
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

void detected(Mat& img, const array<double, 5> & data, int size, 
        const Scalar& color, const double scale=1, int thickness=1, int lineType=8, int shift=0)
{
    for (int direction = 0; direction < 2; direction++)
    {
    const double s = sin(data[2 + direction]);
    const double c = cos(data[2 + direction]);
    double u11 = data[0] - s * data[4];
    double u12 = data[0] + c * size - s * data[4];
    
    double v11 = data[1] + c * data[4];
    double v12 = data[1] + s * size + c * data[4];
    
    line(img, Point(u11*scale + 0.5 * scale, v11*scale+ 0.5 * scale),
            Point(u12*scale+ 0.5 * scale, v12*scale+ 0.5 * scale),
            color, thickness, lineType, shift);
            
    double u21 = data[0] + s * data[4];
    double u22 = data[0] - c * size + s * data[4];
    
    double v21 = data[1] - c * data[4];
    double v22 = data[1] - s * size - c * data[4];
    
    line(img, Point(u21*scale+ 0.5 * scale, v21*scale+ 0.5 * scale),
            Point(u22*scale+ 0.5 * scale, v22*scale+ 0.5 * scale),
            color, thickness, lineType, shift);
    }
}

//CornerDetector function definition

void CornerDetector::improveCorners(Vector2dVec & pointVec) const
{
//    Mat8uc3 img(_img.size());
//    cvtColor(_img, img, CV_GRAY2BGR);
//    resize(img, img, Size(0, 0), 4.5, 4.5);
    
    for (auto & x : pointVec)
    {
        array<double, 5> dataArr;
        initPoin(x, dataArr.data());
        
//        detected(img, dataArr, 15, Scalar(255, 200), 4.5);
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
        /////////////
//        detected(img, dataArr, 15, Scalar(0, 200, 255), 4.5);
        
        
        
//        
//        cout << dataArr[4] << endl;
//        imshow("transitions", img);
//    waitKey();
           
        x[0] = dataArr[0];
        x[1] = dataArr[1];
    }
    
}

CornerDetector::CornerDetector(const Mat8u & img, const int initRadius) : //: _img(img)
    _gradx(img.rows, img.cols),
    _grady(img.rows, img.cols),
    _img(img.rows, img.cols),
    INIT_RADIUS(initRadius)
{
    //TODO instead of derivating the whole image, pass just a small patch to the optimization
    Sobel(img, _gradx, CV_32F, 1, 0, 3, 1e-3);
    Sobel(img, _grady, CV_32F, 0, 1, 3, 1e-3);
    img.copyTo(_img);
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
    Vector2i pt0(pti[0] + INIT_RADIUS, pti[1]);
    CurveRasterizer<int, Polynomial2> raster(pti[0] + INIT_RADIUS, pti[1],
                                             pti[0], pti[1] + INIT_RADIUS, circle);
    //////////
    
    //TODO replace with central differences
    //go along the circle and save the transition strengths
    vector<double> sampleVec;
    vector<double> transitionVec;
    Vector2iVec circleVec;
    for (int i = 0; i < ceil(2 * INIT_RADIUS * M_PI); i++, raster.step())
    {
        if (i > 5) //check whether whe have done a lap
        {
            // the actual point is in the neighborhood of the initial point?
            if (abs(raster.u - pt0[0]) <= 1 and abs(raster.v - pt0[1]) <= 1) break;
        }
        sampleVec.push_back(_img(raster.v, raster.u));
        circleVec.emplace_back(raster.u, raster.v);
    }
    //differentiate using the central differences
    const int sampleCount = sampleVec.size();
    transitionVec.resize(sampleCount);
    for (int i = 0; i < sampleCount; i++)
    {
        int idx1 = (i + sampleCount - 1) % sampleCount;
        int idx2 = (i + 1) % sampleCount;
        transitionVec[i] = 0.5 * (sampleVec[idx2] - sampleVec[idx1]);
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
    data[4] = 0;
}

    
