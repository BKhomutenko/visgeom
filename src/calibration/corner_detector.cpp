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

SubpixelCorner::SubpixelCorner(const Mat32f & gradu, const Mat32f & gradv, const Vector2d pt0,
             const int steps, const double length) :    
            _graduGrid(gradu.cols, gradu.rows, (float*)gradu.data),
            _gradvGrid(gradv.cols, gradv.rows, (float*)gradv.data),
            _stepLength(length / steps),
            _prior(pt0)
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
    
    //prior position
    *cost = 0.1 * ( pow(_prior[0] - u, 2) + pow(_prior[1] - v, 2));
    if (gradient != NULL)
    {
        gradient[0] = 0.2 * (u - _prior[0]);
        gradient[1] = 0.2 * (v - _prior[1]);
        fill(gradient + 2, gradient + 5, 0.);
    }
    
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
    
    line(img, Point(u21*scale + 0.5 * scale, v21*scale + 0.5 * scale),
            Point(u22*scale + 0.5 * scale, v22*scale + 0.5 * scale),
            color, thickness, lineType, shift);
    }
}


//CornerDetector function definition

void CornerDetector::improveCorners(Vector2dVec & pointVec)
{
    vector<double> radVec;
    radVec.reserve(pointVec.size());
    for (int i = 0; i < pointVec.size(); i++)
    {
        //get the maximum radius
        double radMax = 7;
        if (i > _Nx) radMax = min( radMax, (pointVec[i] - pointVec[i - _Nx]).norm()*0.7 );
        else radMax = min( radMax, (pointVec[i] - pointVec[i + _Nx]).norm()*0.7 );
        if (i > 0) radMax = min( radMax, (pointVec[i] - pointVec[i - 1]).norm()*0.7 );
        else radMax = min( radMax, (pointVec[i] - pointVec[i + 1]).norm()*0.7 );
        radVec.push_back(radMax);
    }
    for (int i = 0; i < pointVec.size(); i++)
    {
        Vector2d x = pointVec[i];
        array<double, 5> dataArr;
        initPoin(round(x), dataArr.data());
        ceres::GradientProblem problem(new SubpixelCorner(_gradx, _grady, x, 7, radVec[i]));
        ceres::GradientProblemSolver::Options options;
//        options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
//        options.line_search_type = ceres::ARMIJO;
//        options.line_search_sufficient_function_decrease = 1e-6;
//        options.max_num_line_search_step_size_iterations = 5;
        options.logging_type = ceres::SILENT;
        options.minimizer_progress_to_stdout = true;
        ceres::GradientProblemSolver::Summary summary;
        ceres::Solve(options, problem, dataArr.data(), &summary);
        
//        cout << summary.FullReport() << endl;
        
        pointVec[i][0] = dataArr[0];
        pointVec[i][1] = dataArr[1];
    }
    
}

CornerDetector::CornerDetector(int Nx, int Ny, int initRadius, bool improveDetection, bool debug) :
    _Nx(Nx),
    _Ny(Ny),
    MAX_CANDIDATE_COUNT(10 * Nx * Ny),
    INIT_RADIUS(initRadius),
    IMPROVE_DETECTION(improveDetection),
    DEBUG(debug)
{ }

void CornerDetector::setImage(const Mat8u & img)
{
    img.copyTo(_img);
    
    _resp.create(_img.size());
    _imgrad.create(_img.size());
    _gradx.create(_img.size());
    _grady.create(_img.size());
    _src1.create(_img.size());
    _src2.create(_img.size());
    _detected.create(_img.size());
    _idxMap.create(_img.size());
}

bool CornerDetector::detectPattern(Vector2dVec & ptVec)
{
    const vector<double> sigmaVec = {1.4, 2, 1};
//    const vector<double> sigmaVec = {3};
    vector<int> idxVec;
    bool detected = false;
    for (int i = 0; i < sigmaVec.size() and not detected; i++)
    {
        INIT_RADIUS = round(1.5 * sigmaVec[i]);
        computeResponse(0.7, sigmaVec[i]);
        
        selectCandidates();
        if (DEBUG) cout << "HEAP SIZE : " << _hypHeap.size()  << endl;
        if (_hypHeap.size() < _Nx * _Ny) continue;
        
        constructGraph();
        
        idxVec = selectPattern();
        
        if (idxVec.size() != _Nx * _Ny) continue;
        
        detected = true;
    }
    if (not detected) return false;
    ptVec.clear();
    for (auto & idx : idxVec)
    {
        auto pt = _ptVec[idx];
        ptVec.emplace_back(pt[0], pt[1]);
    }
    
    if (IMPROVE_DETECTION)
    {
        improveCorners(ptVec);
    }
    
    return true;
}
    
void CornerDetector::computeResponse(const double SIGMA_1, const double SIGMA_2)
{
//    const double SIGMA_1 = 0.7; //TODO make a parameter
    const int FILTER_SIZE_1 = 3;
    GaussianBlur(_img, _src1, Size(FILTER_SIZE_1, FILTER_SIZE_1), SIGMA_1, SIGMA_1);
    
//    const double SIGMA_2 = 2.5; //TODO make a parameter
    const int FILTER_SIZE_2 = 1 + 2 * ceil(SIGMA_2);
    GaussianBlur(_img, _src2, Size(FILTER_SIZE_2, FILTER_SIZE_2), SIGMA_2, SIGMA_2);
    _imgrad.setTo(0);
    _resp.setTo(0);
    double acc = 0;
    int count = 0;
    const int SIZE = 1;
    const int HSIZE = 1;
    for (int v = SIZE; v < _img.rows - SIZE; v++)
    {
        for (int u = SIZE; u < _img.cols - SIZE; u++)
        {
            //compute sharp gradient
            double gxSharp = ( _src1(v, u + 1) - _src1(v, u - 1)  
- 0.3*(_src2(v, u + 1) - _src2(v, u - 1)) ) / 2.;
            double gySharp = ( _src1(v + 1, u) - _src1(v - 1, u) 
- 0.3*(_src2(v + 1, u) - _src2(v - 1, u)) ) / 2.;
//            double gxSharp = ( _src1(v, u + 1) - _src1(v, u - 1) ) / 2.;
//            double gySharp = ( _src1(v + 1, u) - _src1(v - 1, u) ) / 2.;
            _gradx(v, u) = gxSharp * 0.01;
            _grady(v, u) = gySharp * 0.01;
            _imgrad(v, u) = sqrt(gxSharp*gxSharp + gySharp*gySharp) * 0.01;
            
            
            
            
            //We are looking for saddle points, that is negative hessian determinant
            double iuu = _src2(v, u - SIZE) + _src2(v, u + SIZE) - 2* _src2(v, u);
            iuu /= SIZE*SIZE;
            double ivv = _src2(v - SIZE, u) + _src2(v + SIZE, u) - 2* _src2(v, u);
            ivv /= SIZE*SIZE;
            double iuv = _src2(v - HSIZE, u - HSIZE) + _src2(v + HSIZE, u + HSIZE) 
                            - _src2(v  + HSIZE, u - HSIZE) - _src2(v - HSIZE, u + HSIZE);
            iuv /= 4*HSIZE*HSIZE;
            
            double gx = (_src2(v, u + HSIZE) - _src2(v, u - HSIZE) ) / (2 * HSIZE);
            double gy = (_src2(v + HSIZE, u) - _src2(v - HSIZE, u) ) / (2 * HSIZE);
            double gsq = gx * gx + gy * gy;
            double _respVal = -iuu * ivv + iuv * iuv - 0.001*pow(gsq, 2);
            if (_respVal > 0.01) 
            {
                _resp(v, u) = _respVal;
                acc += _respVal;
                count++;
            }
        }
    }
    _avgVal = acc / count;
    if (DEBUG)
    {
        imshow("img", _img );
        imshow("imgrad", _imgrad);
        imshow("resp", _resp / 100);
        imwrite("img.png", _img );
        
        Mat32f resp2;
        GaussianBlur(_resp, resp2, Size(11, 11), 1, 1);
        imwrite("resp.png", 255 - 25*resp2);
        waitKey();
    }
}

bool CornerDetector::checkCorner(const Vector2i & pt, const int checkRadius)
{   
//    if (DEBUG) cout << "check corner" << endl;
//    return true;
    int faults = 0;
    const int MAX_FAULTS = 0;
    const int RADIUS_MAX = checkRadius + max(3, checkRadius);
//    cout << checkRadius << "   " << RADIUS_MAX << endl;

    const double ALPHA_1 = 0.3;
    const double ALPHA_2 = 0.5;
    for (int radius = checkRadius; radius < RADIUS_MAX and faults <= MAX_FAULTS; radius++)
    {
//        if (DEBUG) cout << "scale : " << radius << endl;
        Polynomial2 circle = Polynomial2::Circle(pt[0], pt[1], radius);
        
        Vector2i pt0(pt[0] + radius, pt[1]);
        CurveRasterizer<int, Polynomial2> raster(pt[0] + radius, pt[1],
                                                 pt[0], pt[1] + radius, circle);
        
        
        Vector2iVec circleVec = getCircle(pt, radius);
        vector<double> sampleVec = getSamples(circleVec);
        vector<double> transitionVec = centralDifferences(sampleVec);
//       
//        if (DEBUG)
//        { 
//            cout << "faults : " << faults << endl;
//            printVector(transitionVec.begin(), transitionVec.end(), 6);
//            printVector(sampleVec.begin(), sampleVec.end(), 6);
//            Mat8u imgTest;
//            _img.copyTo(imgTest);
//            for (auto & pt : circleVec)
//            {
//                imgTest(pt[1], pt[0]) = 255;
//            }
//            imshow("img", imgTest);
//        }
        
        // find the strongest
        int distThresh = sampleVec.size() / 2 - 2;
        int distThresh2 = sampleVec.size() - distThresh;
        const auto itMax = max_element(transitionVec.begin(), transitionVec.end());
        double transMax = *itMax;
//        cout << "max1 : " << *itMax << endl;
        // remove the maximum and its neighbors
        setZero(transitionVec.begin(), transitionVec.end(), itMax);
        // check the other three : must be at least 0.75
        
        const auto itMax2 = max_element(transitionVec.begin(), transitionVec.end());
//        cout << "max2 : " << *itMax2 << endl;
        
        if (*itMax2 < transMax * ALPHA_1)  { faults++; continue; } //return false;
        double transMax2 = *itMax2;
        setZero(transitionVec.begin(), transitionVec.end(), itMax2);
        
//        cout << "dist : " <<  abs(distance(itMax, itMax2)) << " " << distThresh << endl;
        if (abs(distance(itMax, itMax2)) < distThresh or 
            abs(distance(itMax, itMax2)) > distThresh2) { faults++; continue; } //return false;
            
        auto itMax3 = max_element(transitionVec.begin(), transitionVec.end());
//        cout << "max3 : " << *itMax3 << endl;
        if (*itMax3 > transMax2 * ALPHA_2) { faults++; continue; } //return false;
        setZero(transitionVec.begin(), transitionVec.end(), itMax3);
        
        const auto itMin = min_element(transitionVec.begin(), transitionVec.end());
        double transMin = *itMin;
//        cout << "min1 : " << *itMin << endl;
        
        if (transMin > transMax * -ALPHA_1) { faults++; continue; } //return false;
        setZero(transitionVec.begin(), transitionVec.end(), itMin);
        
        const auto itMin2 = min_element(transitionVec.begin(), transitionVec.end());
//        cout << "min2 : " << *itMin2 << endl;
        
        if (*itMin2 > transMin * -ALPHA_1) { faults++; continue; } //return false;
        double transMin2 = *itMin2;
        setZero(transitionVec.begin(), transitionVec.end(), itMin2);
        
//        cout << abs(distance(itMin, itMin2)) << " " << distThresh << endl;
        if (abs(distance(itMin, itMin2)) < distThresh or
            abs(distance(itMin, itMin2)) > distThresh2) { faults++; continue; } //return false;
        
        auto itMin3 = min_element(transitionVec.begin(), transitionVec.end());
//        cout << "min3 : " << *itMin3 << endl;
        if (*itMin3 < transMin2 * ALPHA_2) { faults++; continue; } //return false;
        setZero(transitionVec.begin(), transitionVec.end(), itMin3);
    }
    if (faults > MAX_FAULTS) return false;
    
    return true;
//    if ( abs( abs(distance(itMin, itMin2)) - abs(distance(itMax, itMax2)) ) > 2 )
//    {
//        if (DEBUG)
//        {
//            for (auto & it : {itMin, itMin2, itMax, itMax2})
//            {
//                int d = distance(transitionVec.begin(), it);
//                _detected(circleVec[d][1], circleVec[d][0]) = 128;
//            }
////            cout << "REJECTED : " << pt.transpose() << endl;
////            cout << abs(distance(itMax, itMax2)) << "   " << abs(distance(itMin, itMin2)) << endl;
//        }
//        return false;
//     
//    }
    
    
    
}

//TODO make the window size adaptive
//and more generally try to detect the board on different scales with different filter size etc
bool CornerDetector::scaleInvarient(const Vector2i & pt)
{
    for (int radius = INIT_RADIUS; radius < INIT_RADIUS * 2; radius++)
    {
        const double MIN_GRAD_THRESH = 1e-3;
        const double CRITERION_THRESH = 0.3;
        double acc = 0;
        int count = 0;
        double normAcc = 1e-10;
        for (int dv = -radius; dv <= radius; dv++)
        {
            for (int du = -radius; du <= radius; du++)
            {
                double sqNorm = du*du + dv*dv;
                if (sqNorm > radius*radius + 1 or sqNorm < 1) continue;
                int u = pt[0] + du;
                int v = pt[1] + dv;
                if (u < 0 or u >= _img.cols or v < 0 or v >= _img.rows) continue;
                double gx = _gradx(v, u);
                double gy = _grady(v, u);
                double gradSqNorm = gx * gx + gy * gy;
                if (gradSqNorm < MIN_GRAD_THRESH) continue;
    //            acc += pow(gx * du + gy * dv, 2) / (sqNorm * gradSqNorm);
    //            count++;
                acc += pow(gx * du + gy * dv, 2) / sqNorm;
                normAcc += gradSqNorm;

            }
            
        }
//        if (DEBUG)
//        {
//            cout << pt.transpose() << endl;
//            cout << "invariance criterion : "  << acc / normAcc << endl;
//            if (acc / normAcc < CRITERION_THRESH)
//            {
//                _detected(pt[1], pt[0]) = 128;
//            }
//            else
//            {
//                _detected(pt[1], pt[0]) = 25;
//            }
//            imshow("detected", _detected);
//            waitKey();
//        }
        if (acc / normAcc < CRITERION_THRESH) return true;
    }
    return false;
}

void CornerDetector::selectCandidates()
{
    //find local maxima
    const int WINDOW_SIZE = INIT_RADIUS;
    vector<pair<double, Vector2i>> maximaHeap;
    for (int v = WINDOW_SIZE; v < _img.rows - WINDOW_SIZE; v++)
    {
        for (int u = WINDOW_SIZE; u < _img.cols - WINDOW_SIZE; u++)
        {
            
            bool isMax = true;
            const double val = _resp(v, u);
            if (val < _avgVal) continue;
            for (int j = -WINDOW_SIZE; j <= WINDOW_SIZE; j++)
            {
                for (int i = -WINDOW_SIZE; i <= WINDOW_SIZE; i++)
                {
                    
                    if (i == 0 and j == 0) continue;
                    double sqNorm = i*i + j*j;
                    if (sqNorm > WINDOW_SIZE*WINDOW_SIZE + 1) continue;
                    const double valNeigh = _resp(v + j, u + i);
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
                maximaHeap.emplace_back(val, Vector2i(u, v));
            }
        }
    }
    
    make_heap(maximaHeap.begin(), maximaHeap.end(), comp);
    
    ///compute thresh
    vector<pair<double, Vector2i>> threshHeap = maximaHeap;
    
    double acc = 0;
    const int REF_ITEM_COUNT = _Nx * _Ny;
    for (int i = 0; i < REF_ITEM_COUNT and not threshHeap.empty(); i++)
    {
        pop_heap(threshHeap.begin(), threshHeap.end(), comp);
        acc += threshHeap.back().first;
        threshHeap.pop_back();
    }
    
    const double VAL_THRESH = 0.05 * acc / REF_ITEM_COUNT;
    
    
    /// Check coner-ness   
    if (DEBUG) _detected.setTo(0);
    
    _hypHeap.clear();   
    int i = 0;                                      
    while (not maximaHeap.empty() and 
            maximaHeap.front().first > VAL_THRESH and 
            i < MAX_CANDIDATE_COUNT)
    {
        pop_heap(maximaHeap.begin(), maximaHeap.end(), comp);
        auto pt = maximaHeap.back().second;
        maximaHeap.pop_back();
        // compute transitions
        
        bool checked = false;
//        if (DEBUG) cout << "START CHECKING" << endl;
        for (int radius = 1; radius <= INIT_RADIUS and not checked; radius++)
        {
            checked = checkCorner(pt, radius);
        }
        
        if (not checked)
        {
            if (DEBUG) 
            {
//                _detected(pt[1], pt[0]) = 64;
                cross(_detected, pt[0], pt[1], 2, 255,3);
//                imshow("detected", _detected);
//                waitKey();
            }
            continue;
        }
        else if (DEBUG) 
        {
//                _detected(pt[1], pt[0]) = 64;
            cross(_detected, pt[0], pt[1], 2, 200, 3);
//                imshow("detected", _detected);
//                waitKey();
        }
        
        if (not scaleInvarient(pt)) continue;
        if (DEBUG)
        {
            cross(_detected, pt[0], pt[1], 7, 150, 3);
        }
        
        _hypHeap.emplace_back(-pt[0] - pt[1], pt);
        i++;
    }
    if (DEBUG)
    {
        imshow("detected", _detected);
        imwrite("detected.png", 255 -  _detected);
        waitKey();
    }
    make_heap(_hypHeap.begin(), _hypHeap.end(), comp);
}


void CornerDetector::constructGraph()
{
    _idxMap.setTo(-1);
    Mat8u timeMat(_idxMap.size());
    timeMat.setTo(0);
    queue<TimePoint> fringe;
    //init the fringe
    //TODO make it after the previous filtering
    _arcVec.clear();
    _ptVec.clear();
    _gradThresh.clear();
    if (DEBUG) _detected.setTo(0);
    for (int i = 0; not _hypHeap.empty(); i++)
    {
        pop_heap(_hypHeap.begin(), _hypHeap.end(), comp);
        Vector2i pt = _hypHeap.back().second;
        _hypHeap.pop_back();
        Vector2iVec initVec = getTransitions(pt);
        _gradThresh.push_back(DOUBLE_MAX);
        for (auto & ptInit : initVec)
        {
            //TODO check the gradient angle 
            fringe.emplace(0, i, ptInit[0], ptInit[1]);
            double gradAbs = _imgrad(ptInit[1], ptInit[0]);
            _gradThresh.back() = min(gradAbs / 2 , _gradThresh.back());
            
        }
        if (DEBUG)
        {
            for (auto & ptt : initVec)
            {
                _detected(ptt[1], ptt[0]) = 180;
                cout << " grad : " << _imgrad(ptt[1], ptt[0])  << "   " << _gradThresh.back() << endl;
                cout << ptt << endl;
            }
//            imshow("init", _detected);
//            waitKey();
        }
        
        _arcVec.emplace_back();
        _ptVec.push_back(pt);
    }
    
    if (DEBUG)
    {
        imshow("init", _detected);
        waitKey();
    }
        
    const vector<int> duVec = {-1, 0, 1, 1, 1, 0, -1, -1};
    const vector<int> dvVec = {-1, -1, -1, 0, 1, 1, 1, 0};
    
    const int CANDIDATE_COUNT = _arcVec.size();
    const int SEARCH_REACH = 140;
    const double SQUARED_GRAD_MAX = 5;
    //TODO make adaptive
    int iter = 0;
    while (not fringe.empty() and fringe.front().t < SEARCH_REACH)
    {
        iter++;
        auto e = fringe.front();
        fringe.pop();
        if (_idxMap(e.v, e.u) != -1)
        {
            continue;            
        }
        _idxMap(e.v, e.u) = e.idx;
        if (DEBUG) timeMat(e.v, e.u) = e.t + 150;
        //first connect
        bool connected = false;
        for (int i = 0; i < 8; i++)
        {
            const int u2 = e.u + duVec[i];
            const int v2 = e.v + dvVec[i];
            if (u2 < 0 or u2 >= _idxMap.cols or v2 < 0 or v2 >= _idxMap.rows)
            {
                continue;
            }
            const int idx2 = _idxMap(v2, u2);
            if (idx2 != -1 and idx2 != e.idx)
            {
                connected = true;
                if (find(_arcVec[e.idx].begin(), _arcVec[e.idx].end(), idx2) == _arcVec[e.idx].end())
                {
                    //Computing the arc sign
                    Vector2d arc = (_ptVec[idx2] - _ptVec[e.idx]).cast<double>();
                    Vector2d normalArc = arc.normalized();
                    int signAcc = 0;
                    for (int base = 1; base <= INIT_RADIUS; base++)
                    {
//                        cout << base << endl;
                        Vector2d stepVec = normalArc * base;
                        const int LAMBDA_MAX = 5;
                        for (int lambda = 1; lambda < LAMBDA_MAX; lambda++)
                        {
                            Vector2d midPoint = _ptVec[e.idx].cast<double>() + arc * (double(lambda) / LAMBDA_MAX);
                            double sample1 = bilinear<double>(_src2, midPoint[0] - stepVec[1],
                                                                     midPoint[1] + stepVec[0]);
                            double sample2 = bilinear<double>(_src2, midPoint[0] + stepVec[1], 
                                                                     midPoint[1] - stepVec[0]);
                            signAcc += sign(sample1 - sample2);
                        }
                    }
                    signAcc = sign(signAcc);
            
                    
                    _arcVec[e.idx].push_back(idx2);
                    _arcVec[idx2].push_back(e.idx);
                    _arcSign[make_pair(e.idx, idx2)] = signAcc;
                    _arcSign[make_pair(idx2, e.idx)] = -signAcc;
                }
            }
        }    
            
        if (not connected)
        {
            double bestGradProj = 0;
            int bestDir;
            for (int i = 0; i < 8; i++)
            {
                //proliferate
                const int u2 = e.u + duVec[i];
                const int v2 = e.v + dvVec[i];
                if (u2 < 0 or u2 >= _idxMap.cols or v2 < 0 or v2 >= _idxMap.rows)
                {
                    continue;
                }
                if (_idxMap(v2, u2) != -1) continue;
                
                double x1 = u2 - _ptVec[e.idx][0];
                double y1 = v2 - _ptVec[e.idx][1];
                double t = sqrt(x1 * x1 + y1 * y1);
    //            if (e.u == 438 and e.v == 129) cout << "imgrad : " << _imgrad(v2, u2) << endl;
                double gradProj = abs(_gradx(v2, u2) * y1 - _grady(v2, u2) * x1);
                
                if ( t > 0 and gradProj / t  < _gradThresh[e.idx])
                {
                     continue;
                }
    //            else if (DEBUG)
    //            {
    //                cout << t2 << " " << _imgrad(v2, u2) << " " << _gradThresh[e.idx] << endl;
    //            }
                if (t > 0)
                {
                    
    //                check the gradient orientation
                    double xg = _gradx(v2, u2);
                    double yg = _grady(v2, u2);
                    double c = x1 * xg + y1 * yg;
                    double s = x1 * yg - y1 * xg;
                    double angle = abs(atan2(s, c));
                    double thresh = min(M_PI / 5, t / 150. + 1 / t);
    //                if (e.u == 438 and e.v == 129) cout << M_PI / 2 - thresh << "   " 
    //                                                    << angle << "   "
    //                                                    << M_PI / 2 + thresh << "   " << endl;
                    
                    if (angle < M_PI / 2 - thresh or 
                        angle > M_PI / 2 + thresh) continue;
                        
//                    if (angle < M_PI / 3 or 
//                        angle > M_PI * 2 / 3 ) continue;
                }
                
                fringe.emplace(e.t + 1, e.idx, u2, v2);
            }
        }
//        if (DEBUG and iter % 100 == 0)
//        {
//            Mat16s idxView;
//            resize(_idxMap, idxView, Size(0, 0), 4, 4, cv::INTER_NEAREST);
//            for (auto & x : idxView)
//            {
//                x = (x * 2000) % 32000;
//            }    
//            
//            imshow("detected", idxView  - 32000);
//            waitKey();
//        }
    }
    
    if (DEBUG)
    {
        Mat8u graphView;
        int RESIZE_SCALE = 1;
        resize(_detected, graphView, Size(0, 0), RESIZE_SCALE, RESIZE_SCALE, cv::INTER_NEAREST);   
        graphView.setTo(0);
        for (auto & pt : _ptVec)
        {
            cross(graphView, pt[0], pt[1], 7, 180, 3);
        }
        for (int i = 0; i < _arcVec.size(); i++)
        {
            for (auto arc : _arcVec[i])
            {
                if (arc < i) continue;
                int x1 = _ptVec[i][0] * RESIZE_SCALE;
                int y1 = _ptVec[i][1] * RESIZE_SCALE;
                int x2 = _ptVec[arc][0] * RESIZE_SCALE;
                int y2 = _ptVec[arc][1] * RESIZE_SCALE;
                
                line(graphView, Point(x1, y1),
                    Point(x2, y2),
                    255, 2, 8, 0);
                    
//                if (_arcSign[make_pair(i, arc)] == 1)
//                {
//                    cross(graphView, (x1 + 3 * x2) / 4,(y1 + 3 * y2) / 4, 2, 55, 3);
//                }
//                else  cross(graphView, (x2 + 3 * x1) / 4,(y2 + 3 * y1) / 4, 2, 55, 3);
            }
        }
        
        
        
        imshow("graph", graphView);
        imshow("timeMat", timeMat); 
        
        imwrite("graph.png", 255 -  graphView);
        imwrite("timeMat.png", 255 - timeMat); 
        
        waitKey();
    }
    
}

double compareVectors(Vector2i v1, Vector2i v2)
{
    return (v1 - v2).norm() / double(v1.norm());
//    double norm1 = v1.norm();
//    double norm2 = v2.norm();
//    double dlen = (norm2 - norm1) / norm1;
//    double s = v1[0]*v2[1] - v1[1]*v2[0];
//    double c = v1.dot(v2);
//    double th = atan2(s, c);
//    return abs(th); //FIXME parameters to tune
}

vector<int> CornerDetector::extractSequence(int idx0, int idx1)
{
    vector<int> chain;
    
    if (DEBUG) 
    {
        cout << "NEW CHAIN EXTRACTION" << endl;
        cout << _ptVec[idx0].transpose() << "    " << _ptVec[idx1].transpose() << endl;
    }
    ///FIND A GOOD DIRECTIOIN///
    
    Vector2i d1 = _ptVec[idx1] - _ptVec[idx0];
    
    int idx2 = -1;
    double bestNormDiff = 1; //TODO a parameter to be set
//    cout << _arcVec[idx1].size() << endl;
    for (auto & n2 : _arcVec[idx1])
    {
//        cout << 
        
        
        Vector2i d2 = _ptVec[n2] - _ptVec[idx1];
        double normDiff = compareVectors(d1, d2);
        if (DEBUG) cout << d1.transpose() << "   " << d2.transpose() << endl << normDiff << endl;
        if (_arcSign[make_pair(idx1, n2)] == _arcSign[make_pair(idx0, idx1)]) continue;
        if (normDiff < bestNormDiff)
        {
            bestNormDiff = normDiff;
            idx2 = n2;
        }
    }
    if (idx2 == -1) return chain;
    
    chain.push_back(idx0);
    chain.push_back(idx1);
    chain.push_back(idx2);

    ///EXTRACT THE CHAIN///
    while (true)
    {   
        if (DEBUG) cout << "New step" << endl;
        auto chainIt = chain.end();
        int idx1 = *(--chainIt);
        int idx0 = *(--chainIt);
        Vector2i d0 = _ptVec[idx1] - _ptVec[idx0];
        if (DEBUG) cout << d0.transpose() << endl;
        int idx2 = -1;
        double bestNormDiff = 1;
        for (auto & n2 : _arcVec[idx1])
        {
            Vector2i d1 = _ptVec[n2] - _ptVec[idx1];
            if (DEBUG) cout << d1.transpose() << endl;
            if (n2 == idx0) continue;
            if (_arcSign[make_pair(idx1, n2)] == _arcSign[make_pair(idx0, idx1)]) continue;
            
            double normDiff = compareVectors(d0, d1);
     if (DEBUG)       cout << normDiff << endl;
            if (normDiff < bestNormDiff)
            {
                bestNormDiff = normDiff;
                idx2 = n2;
            }
        } 
//        cout << endl;
        if (idx2 != -1)
        {
            chain.push_back(idx2);
        }
        else break;
    }
    
    return chain;
}

//TODO implement
//vector<int> CornerDetector::weaveGrid(const vector<int> & chainX, const vector<int> & chainY)
//{

//}

vector<int> CornerDetector::selectBestOrthogonalChain(const int idx0, const int idx1,
                                                     const int EPS, const int LENGTH)
{
    double bestCost = 0.3;
    vector<int> bestChain;
    int baseSign = _arcSign[make_pair(idx0, idx1)];
    for (auto & nx : _arcVec[idx0])
    {
        if (nx == idx1) continue;
        if (_arcSign[make_pair(idx0, nx)] == baseSign) continue;
        //check the orientation
        Vector2i d1 = _ptVec[idx1] - _ptVec[idx0];
        Vector2i d2 = _ptVec[nx] - _ptVec[idx0];
        
//        if (DEBUG)
//        {
//            cross(chainMat, _ptVec[nx][0], _ptVec[nx][1], 3, 255, 3);
//            imshow("chainMat", chainMat);
//            waitKey();
//        }
//        cout << _ptVec[idxCur].transpose() << "     " << _ptVec[nx].transpose()  << endl;
        vector<int> chain = extractSequence(idx0, nx);
        double cost = EPS*(d1[0] * d2[1] - d1[1] * d2[0]) / double(d1.norm() * d2.norm());
        if (cost < bestCost) continue;
//        if (DEBUG)
//        {
//            for (auto & pt : chain)
//            {
//                cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255);
//            }
//            imshow("chainMat", chainMat);
//            waitKey();
//        }
            
        if (chain.size() >= LENGTH)
        {
            bestCost = cost;
            bestChain = chain;
        }
    }
    while(bestChain.size() > LENGTH) bestChain.pop_back();
    return bestChain;
}

//TODO rewrite with a finite automaton
vector<int> CornerDetector::selectPattern()
{
    vector<int> res;
    for (int idx0 = 0; idx0 < _ptVec.size(); idx0++)
    {
        //FIXME for debug
        Mat8u chainMat;
        if (DEBUG)
        {
            chainMat.create(_detected.size());
            chainMat.setTo(0);
            for (auto & pt : _ptVec)
            {
                cross(chainMat, pt[0], pt[1], 4, 100, 4);
            }
            cross(chainMat, _ptVec[idx0][0], _ptVec[idx0][1], 5, 255, 5);
//            for (auto & n : _arcVec[idx0]) cout << n << "   ";
//            cout << endl;
        }
             
        if (_arcVec[idx0].size() < 2) continue;
        //first, look for a potential upper left corner
        vector<int> chainY;
        vector<int> chainX;
        for (auto & n : _arcVec[idx0])
        {
            vector<int> chain = extractSequence(idx0, n);
            
            if (chain.size() >= _Ny) 
            {
                while (chain.size() > _Ny) chain.pop_back();
                chainX = selectBestOrthogonalChain(idx0, n, -1, _Nx);
                if (chainX.size() == _Nx)
                {
                    chainY = chain;
                    break;
                }
            }
        }
        if (chainY.empty() or chainX.empty()) continue;
        
        if (DEBUG)
        {
            cout << "Chain1 detected" << endl;
            for (auto & pt : chainX)
            {
                cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255, 2);
            }
            for (auto & pt : chainY)
            {
                cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255, 2);
            }
            imshow("chainMat", chainMat);
            waitKey();
            
        }
        
        res = chainX;
        
        for (int i = 1; i < chainY.size(); i++)
        {
            int idxPrev = chainY[i-1];
            int idxCur = chainY[i];
            int prevSign = _arcSign[make_pair(idxCur, idxPrev)];
            double bestCost = 0.3;
            vector<int> bestChain = selectBestOrthogonalChain(idxCur, idxPrev, 1, _Nx);
            if (bestChain.size() >= _Nx)
            {
                res.insert(res.end(), bestChain.begin(), bestChain.begin() + _Nx);
            }
            else break;
        }
        if (res.size() == _Nx * _Ny and verifyDetection(res))
        {
            return res;
        }
        else
        {
            res.clear();
        }
    }    
    
    return res;
}

bool CornerDetector::verifyDetection(const vector<int> & idxVec)
{
    if (idxVec.size() != _Nx * _Ny) return false;
    for (int i = 1; i < _Nx; i++)
    {
        auto begin = _arcVec[idxVec[i]].begin();
        auto end = _arcVec[idxVec[i]].end();
        if (find(begin, end, idxVec[_Nx + i]) == end) return false;
        vector<int> chain = extractSequence(idxVec[i], idxVec[_Nx + i]);
        if (chain.size() < _Ny) return false;
        for (int j = 2; j < _Ny; j++)
        {
            if (chain[j] != idxVec[j*_Nx + i]) return false;
        }
    }
    return true;
}


Vector2iVec CornerDetector::getCircle(const Vector2i & pt, const int radius) const
{
    Vector2iVec res;
    if (radius == 1)
    {
        const vector<int> duVec = {1, 1, 0, -1, -1, -1, 0, 1};
        const vector<int> dvVec = {0, 1, 1, 1, 0, -1, -1, -1};
        for (int i = 0; i < 8; i++)
        {
            res.emplace_back(normalizePoint(Vector2i(pt[0] + duVec[i], pt[1] + dvVec[i])));
        }
    }
    else
    {
        Polynomial2 circle = Polynomial2::Circle(pt[0], pt[1], radius);
    
        Vector2i pt0(pt[0] + radius, pt[1]);
        CurveRasterizer<int, Polynomial2> raster(pt[0] + radius, pt[1],
                                             pt[0], pt[1] + radius, circle);
        int i = 0;
        while (true)
        {
            res.emplace_back( normalizePoint(Vector2i(raster.u, raster.v)) );
            if (i > 5 and abs(raster.u - pt0[0]) <= 1 and abs(raster.v - pt0[1]) <= 1) break;
            i++; 
            raster.step();
        }  
    }
    return res;
}

vector<double> CornerDetector::getSamples(const Vector2iVec & ptVec) const
{
    vector<double> res;
    res.reserve(ptVec.size());
    for (auto & pt : ptVec)
    {
        res.push_back(_img(pt[1], pt[0]));
    }
    return res;
}

vector<double> CornerDetector::centralDifferences(const vector<double> & sampleVec) const
{
    vector<double> res;
    res.reserve(sampleVec.size());
    res.push_back(sampleVec[1] - sampleVec.back());
    for (int i = 1; i < sampleVec.size() - 1; i++)
    {
        res.push_back(sampleVec[i + 1] - sampleVec[i - 1]);
    }
    int idxLast = sampleVec.size() - 2;
    res.push_back(sampleVec.front() - sampleVec[idxLast]);
    return res;
}

Vector2iVec CornerDetector::getTransitions(const Vector2i & pt)
{
    //find the threshold
//    double gradThresh = 0;
//    for (int j = -INIT_RADIUS; j <= INIT_RADIUS; j++)
//    {
//        for (int i = -INIT_RADIUS; i <= INIT_RADIUS; i++)
//        {
//            
//            if (i == 0 and j == 0) continue;
//            double sqNorm = i*i + j*j;
//            if (sqNorm > INIT_RADIUS*INIT_RADIUS + 1) continue;
//            //FIXME the 100 coef must be replaced in future
//            if (pt[1] + j < 0 or pt[1] + j >= _imgrad.rows or 
//                pt[0] + i < 0 or pt[0] + i >= _imgrad.cols) continue;
//            gradThresh = max(_imgrad(pt[1] + j, pt[0] + i) * 100., gradThresh);
//            
//        }
//    }
//    gradThresh /= 2; //TODO an arbitrary constant
    // for init_radius look for transitions
    Vector2iVec res;
    bool detected = false;
    double bestMaxTransition;
    for (int radius = 1; radius <= INIT_RADIUS + 1; radius++)
    {
        
        //go along the circle and save the transition strengths
        
        Vector2iVec circleVec = getCircle(pt, radius);
        vector<double> sampleVec = getSamples(circleVec);
        vector<double> transitionVec = centralDifferences(sampleVec);
        
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
        if (detected)
        {
            if (bestMaxTransition > 0.7 * transitionVec[maxIdx1]) break;
            else 
            {
                detected = false;
                res.clear();
            }
        }
//        if (transitionVec[maxIdx1] < gradThresh) continue;
        if (transitionVec[maxIdx2] < 0.4 * transitionVec[maxIdx1]) continue;
        else if (transitionVec[minIdx2] > 0.4 * transitionVec[minIdx1]) continue;
        else
        {
            //compute the angles and the intersection of the lines
            res.reserve(4);
            res.push_back(circleVec[maxIdx1]);
            res.push_back(circleVec[maxIdx2]);
            res.push_back(circleVec[minIdx1]);
            res.push_back(circleVec[minIdx2]);
            bestMaxTransition = transitionVec[maxIdx1];
            detected = true;
//            if (DEBUG)
//            {
//                for (auto & ptt : res)
//                {
//                    _detected(ptt[1], ptt[0]) = 180;
//                }
//                imshow("detected", _detected);
//                waitKey();
//            }
        }
    }
    return res;
}

int CornerDetector::initPoin(const Vector2i & pt, double * data)
{
   
    //compute the angles and the intersection of the lines
    Vector2iVec transitions = getTransitions(pt);
    const Vector2i & A = transitions[0];    // max1
    const Vector2i & C = transitions[1];    // max2
    const Vector2i & B = transitions[2];    // min1
    const Vector2i & D = transitions[3];    // min2
    
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
    
    
//    if (DEBUG)
//    {
//        cout << M << endl;
//        cout << b.transpose() << endl;
//        cout << M.inverse() << endl;
//        printPointVector(transitions.begin(), transitions.end());
//        printVector(data, data + 5, 12);
//    }
}

    
