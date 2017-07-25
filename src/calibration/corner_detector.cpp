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
    
    line(img, Point(u21*scale+ 0.5 * scale, v21*scale+ 0.5 * scale),
            Point(u22*scale+ 0.5 * scale, v22*scale+ 0.5 * scale),
            color, thickness, lineType, shift);
    }
}


//CornerDetector function definition

void CornerDetector::improveCorners(Vector2dVec & pointVec) const
{
    
    for (auto & x : pointVec)
    {
        array<double, 5> dataArr;
        initPoin(x, dataArr.data());
        
        ceres::GradientProblem problem(new SubpixelCorner(_gradx, _grady, x));
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
        
        x[0] = dataArr[0];
        x[1] = dataArr[1];
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
    
    computeResponse();
    
    selectCandidates();
    
    if (_hypHeap.size() < _Nx * _Ny) return false;
    
    constructGraph();
    
    vector<int> idxVec = selectPattern();
    
    if (idxVec.size() != _Nx * _Ny) return false;
    
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
    
void CornerDetector::computeResponse()
{
    const double SIGMA_1 = 0.7; //TODO make a parameter
    const int FILTER_SIZE_1 = 3;
    GaussianBlur(_img, _src1, Size(FILTER_SIZE_1, FILTER_SIZE_1), SIGMA_1, SIGMA_1);
    
    const double SIGMA_2 = 3.5; //TODO make a parameter
    const int FILTER_SIZE_2 = 1 + 2 * round(SIGMA_2);
    GaussianBlur(_img, _src2, Size(FILTER_SIZE_2, FILTER_SIZE_2), SIGMA_2, SIGMA_2);
    _imgrad.setTo(0);
    _resp.setTo(0);
    double acc = 0;
    int count = 0;
    const int SIZE = 3;
    const int HSIZE = 2;
    for (int v = SIZE; v < _img.rows - SIZE; v++)
    {
        for (int u = SIZE; u < _img.cols - SIZE; u++)
        {
            //compute sharp gradient
//            double gxSharp = ( _src1(v, u + 1) - _src1(v, u - 1)  - 0.3*(_src2(v, u + 1) - _src2(v, u - 1)) ) / 2.;
//            double gySharp = ( _src1(v + 1, u) - _src1(v - 1, u) - 0.3*(_src2(v + 1, u) - _src2(v - 1, u)) ) / 2.;
            double gxSharp = ( _src1(v, u + 1) - _src1(v, u - 1) ) / 2.;
            double gySharp = ( _src1(v + 1, u) - _src1(v - 1, u) ) / 2.;
            _gradx(v, u) = gxSharp * 0.01;
            _grady(v, u) = gySharp * 0.01;
            _imgrad(v, u) = gxSharp*gxSharp + gySharp*gySharp;
            
            
            
            
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
        imshow ("img", _img );
        imshow ("imgrad", _imgrad / 500);
        imshow ("resp", _resp / 100);
        waitKey();
    }
}

bool CornerDetector::checkCorner(const Vector2i & pt)
{   
//    return true;
    Polynomial2 circle = Polynomial2::Circle(pt[0], pt[1], INIT_RADIUS);
    
    Vector2i pt0(pt[0] + INIT_RADIUS, pt[1]);
    CurveRasterizer<int, Polynomial2> raster(pt[0] + INIT_RADIUS, pt[1],
                                             pt[0], pt[1] + INIT_RADIUS, circle);
    
    vector<double> sampleVec;
    Vector2iVec circleVec;
    for (int i = 0;
         not (i > 5 and abs(raster.u - pt0[0]) <= 1 and abs(raster.v - pt0[1]) <= 1);
         i++, raster.step())
    {
        sampleVec.push_back(_img(raster.v, raster.u));
        circleVec.emplace_back(raster.u, raster.v);
    }
    
    vector<double> transitionVec;
    
    transitionVec.push_back(sampleVec[1] - sampleVec.back());
    
    for (int i = 1; i < sampleVec.size() - 1; i++)
    {
        transitionVec.push_back(sampleVec[i + 1] - sampleVec[i - 1]);
    }
    int idxLast = sampleVec.size() - 2;
    transitionVec.push_back(sampleVec.front() - sampleVec[idxLast]);
                            
    // find the strongest
    int distThresh = sampleVec.size() / 2 - 3;
    
    const auto itMax = max_element(transitionVec.begin(), transitionVec.end());
    double transMax = *itMax;
//        cout << endl << setw(10) << *itMax;
    // remove the maximum and its neighbors
    setZero(transitionVec.begin(), transitionVec.end(), itMax);
    // check the other three : must be at least 0.75
    
    const auto itMax2 = max_element(transitionVec.begin(), transitionVec.end());
//        cout << setw(10) << *itMax2;
    
    if (*itMax2 < transMax * 0.5) return false;
    setZero(transitionVec.begin(), transitionVec.end(), itMax2);
    
    if (abs(distance(itMax, itMax2)) < distThresh) return false;
        
    auto itMax3 = max_element(transitionVec.begin(), transitionVec.end());
    
    if (*itMax3 > transMax * 0.5) return false;
    setZero(transitionVec.begin(), transitionVec.end(), itMax3);
    
    const auto itMin = min_element(transitionVec.begin(), transitionVec.end());
//        cout << setw(10) << *itMin;
    
    if (*itMin > transMax * -0.5) return false;
    setZero(transitionVec.begin(), transitionVec.end(), itMin);
    
    const auto itMin2 = min_element(transitionVec.begin(), transitionVec.end());
//        cout << setw(10) << *itMin2;
    
    if (*itMin2 > transMax * -0.5) return false;
    setZero(transitionVec.begin(), transitionVec.end(), itMin2);
    
    if ( abs(distance(itMin, itMin2)) < distThresh) return false;
    
    auto itMin3 = min_element(transitionVec.begin(), transitionVec.end());
    
    if (*itMin3 < transMax * -0.5) return false;
    setZero(transitionVec.begin(), transitionVec.end(), itMin3);
    
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
    
    
    return scaleInvarient(pt);
}

bool CornerDetector::scaleInvarient(const Vector2i & pt)
{
    const double MIN_GRAD_THRESH = 15e-4;
    const double CRITERION_THRESH = 0.25;
    const int WINDOW_SIZE = 6;
    double acc = 0;
    int count = 0;
    double normAcc = 0;
    for (int dv = -WINDOW_SIZE; dv <= WINDOW_SIZE; dv++)
    {
        for (int du = -WINDOW_SIZE; du <= WINDOW_SIZE; du++)
        {
            double sqNorm = du*du + dv*dv;
            if (sqNorm > WINDOW_SIZE*WINDOW_SIZE + 1 or sqNorm < 5) continue;
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
    if (DEBUG)
    {
        cout << "invariance criterion : "  << acc / normAcc << endl;
    }
    return (acc / normAcc < CRITERION_THRESH);
}

void CornerDetector::selectCandidates()
{
    //find local maxima
    const int WINDOW_SIZE = 5;
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
    
    const double VAL_THRESH = 0.01 * acc / REF_ITEM_COUNT;
    
    
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
        
        if (not checkCorner(pt)) continue;
        
        if (DEBUG) _detected(pt[1], pt[0]) = 255;
        
        _hypHeap.emplace_back(-pt[0] - pt[1], pt);
        i++;
    }
    if (DEBUG)
    {
        imshow ("detected", _detected);
        waitKey();
    }
    make_heap(_hypHeap.begin(), _hypHeap.end(), comp);
}


void CornerDetector::constructGraph()
{
    _idxMap.setTo(-1);

    queue<TimePoint> fringe;
    //init the fringe
    //TODO make it after the previous filtering
    _arcVec.clear();
    _ptVec.clear();
    for (int i = 0; not _hypHeap.empty(); i++)
    {
        pop_heap(_hypHeap.begin(), _hypHeap.end(), comp);
        Vector2i pt = _hypHeap.back().second;
        _hypHeap.pop_back();
        const int SPREAD = 2;
        for (int du = -SPREAD; du <= SPREAD; du += SPREAD)
        {
            for (int dv = -SPREAD; dv <= SPREAD; dv += SPREAD)
            {
                fringe.emplace(0, i, pt[0] + du, pt[1] + dv);
            }
        }
        fringe.emplace(0, i, pt[0], pt[1]);
        _arcVec.emplace_back();
        _ptVec.push_back(pt);
    }
    
    const vector<int> duVec = {-1, 0, 1, 1, 1, 0, -1, -1};
    const vector<int> dvVec = {-1, -1, -1, 0, 1, 1, 1, 0};
    
    const int CANDIDATE_COUNT = _arcVec.size();
    const int SEARCH_REACH = 140;
    const double SQUARED_GRAD_MAX = 15;
    while (not fringe.empty() and fringe.front().t < SEARCH_REACH)
    {
        auto e = fringe.front();
        fringe.pop();
        if (_idxMap(e.v, e.u) != -1)
        {
            continue;            
        }
        _idxMap(e.v, e.u) = e.idx;
        for (int i = 0; i < 8; i++)
        {
            const int u2 = e.u + duVec[i];
            const int v2 = e.v + dvVec[i];
            if (u2 < 0 or u2 >= _idxMap.cols or v2 < 0 or v2 >= _idxMap.rows)
            {
                continue;
            }
            
            
            double x1 = u2 - _ptVec[e.idx][0];
            double y1 = v2 - _ptVec[e.idx][1];
            double t2 = x1 * x1 + y1 * y1;
            if (t2 > 50)
            {
                if (_imgrad(v2, u2) < SQUARED_GRAD_MAX ) continue;
                //check the gradient orientation
                double xg = _gradx(v2, u2);
                double yg = _grady(v2, u2);
                double c = x1 * xg + y1 * yg;
                double s = x1 * yg - y1 * xg;
                double angle = abs(atan2(s, c));
                double t = sqrt(t2);
                double thresh = min(M_PI / 8, t / 150.) + 1. / t;
                if (angle < M_PI / 2 - thresh or 
                    angle > M_PI / 2 + thresh) continue;
            }
            if (_idxMap(v2, u2) == -1)
            {
                fringe.emplace(e.t + 1, e.idx, u2, v2);
            }
            else if (_idxMap(v2, u2) != e.idx)
            {
                const int idx2 = _idxMap(v2, u2);
                if (find(_arcVec[e.idx].begin(), _arcVec[e.idx].end(), idx2) == _arcVec[e.idx].end())
                {
                    _arcVec[e.idx].push_back(idx2);
                    _arcVec[idx2].push_back(e.idx);
                }
            }
        }
    }
    if (DEBUG)
    {
        for (auto & x : _idxMap)
        {
            x = (x * 2000) % 32000;
        }    
        imshow ("detected", _idxMap  - 32000);
    }
}

double compareVectors(Vector2i v1, Vector2i v2)
{
    double norm1 = v1.norm();
    double norm2 = v2.norm();
    double len2 = pow((norm1 - norm2) / norm1, 2);
    double c = v1.dot(v2) / norm1 / norm2;
    double th2 = 2. * (1. - c);
    return th2 + len2/10; //FIXME parameters to tune
}

vector<int> CornerDetector::extractSequence(int idx0, int idx1)
{
    vector<int> chain;
    
    ///FIND A GOOD DIRECTIOIN///
    
    Vector2i d1 = _ptVec[idx1] - _ptVec[idx0];
    
    int idx2 = -1;
    double bestNormDiff = 0.1; //TODO a parameter to be set
    for (auto & n2 : _arcVec[idx1])
    {
        Vector2i d2 = _ptVec[n2] - _ptVec[idx1];
        double normDiff = compareVectors(d1, d2);
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
        auto chainIt = chain.end();
        int idx2 = *(--chainIt);
        int idx1 = *(--chainIt);
        int idx0 = *(--chainIt);
        Vector2i d0 = _ptVec[idx1] - _ptVec[idx0];
        Vector2i d1 = _ptVec[idx2] - _ptVec[idx1];
        Vector2i predicted = 2 * d1 - d0;
        int idx3 = -1;
        double bestNormDiff = 0.1;
//        cout << " predictive search : ";
        for (auto & n2 : _arcVec[idx2])
        {
            Vector2i d2 = _ptVec[n2] - _ptVec[idx2];
            double normDiff = compareVectors(d1, d2);
            if (normDiff < bestNormDiff)
            {
                bestNormDiff = normDiff;
                idx3 = n2;
                break;
            }
        } 
//        cout << endl;
        if (idx3 != -1)
        {
            chain.push_back(idx3);
        }
        else break;
    }
    
    return chain;
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
            _detected.copyTo(chainMat);
            cross(chainMat, _ptVec[idx0][0], _ptVec[idx0][1], 6, 255, 5);
            for (auto & n : _arcVec[idx0]) cout << n << "   ";
            cout << endl;
        }
             
        if (_arcVec[idx0].size() < 2) continue;
        //first, look for a potential upper left corner
        vector<int> chainY;
        vector<int> chainX;
        for (auto & n : _arcVec[idx0])
        {
            vector<int> chain = extractSequence(idx0, n);
            if (DEBUG)
            {
                cout << "Chain1 detected" << endl;
                for (auto & pt : chain)
                {
                    cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255);
                }
                imshow("chainMat", chainMat);
                waitKey();
                
            }
            
            if (chain.size() >= _Ny) 
            {
                while (chain.size() > _Ny) chain.pop_back();
                
                //check whether x-chain is extractable
                for (auto & nx : _arcVec[idx0])
                {
                    if (nx == n) continue;
                    //check the orientation
                    Vector2i d1 = _ptVec[nx] - _ptVec[idx0];
                    Vector2i d2 = _ptVec[n] - _ptVec[idx0];
                    
                    if (d1[0] * d2[1] - d1[1] * d2[0] < 0.1 * d1.norm() * d2.norm()) continue;
                    vector<int> chain2 = extractSequence(idx0, nx);
                    if (chain2.size() >= _Nx)
                    {
                        while (chain2.size() > _Nx) chain2.pop_back();
                        
                        
                        if (DEBUG)
                        {
                            cout << "Chain2 detected" << endl;
                            cout << d1.transpose() << "   " << d2.transpose() << endl;
                            for (auto & pt : chain2)
                            {
                                cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255);
                            }
                            imshow("chainMat", chainMat);
                            waitKey();
                        }
                
                
                
                        chainX = chain2;
                        break;
                    }
                }
                if (chainX.size() == _Nx)
                {
                    chainY = chain;
                    break;
                }
            }
        }
        if (chainY.empty() or chainX.empty()) continue;
        
        res = chainX;
        
        for (int i = 1; i < chainY.size(); i++)
        {
            int idxPrev = chainY[i-1];
            int idxCur = chainY[i];
            bool detected = false;
            for (auto & nx : _arcVec[idxCur])
            {
                if (nx == idxPrev) continue;
                //check the orientation
                Vector2i d1 = _ptVec[idxPrev] - _ptVec[idxCur];
                Vector2i d2 = _ptVec[nx] - _ptVec[idxCur];
                
                if (DEBUG)
                {
                    cross(chainMat, _ptVec[nx][0], _ptVec[nx][1], 3, 255, 3);
                    cout << d1[0] * d2[1] - d1[1] * d2[0] << endl;
                    imshow("chainMat", chainMat);
                    waitKey();
                }
                if (d1[0] * d2[1] - d1[1] * d2[0] < 0.1 * d1.norm() * d2.norm()) continue;
                vector<int> chain2 = extractSequence(idxCur, nx);
                
                if (DEBUG)
                {
                    for (auto & pt : chain2)
                    {
                        cross(chainMat, _ptVec[pt][0], _ptVec[pt][1], 5, 255);
                    }
                    imshow("chainMat", chainMat);
                    waitKey();
                }
                    
                if (chain2.size() >= _Nx)
                {
                    
                    detected = true;
                    res.insert(res.end(), chain2.begin(), chain2.begin() + _Nx);
                    break;
                }
            }
            if (not detected) break;
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
        vector<int> chain = extractSequence(idxVec[i], idxVec[_Nx + i]);
        if (chain.size() < _Ny) return false;
        for (int j = 2; j < _Ny; j++)
        {
            if (chain[j] != idxVec[j*_Nx + i]) return false;
        }
    }
    return true;
}


void CornerDetector::initPoin(const Vector2d & pt, double * data) const
{
    // for init_radius look for transitions
    
    Polynomial2 circle = Polynomial2::Circle(pt[0], pt[1], INIT_RADIUS);
    
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

    
