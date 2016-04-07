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
Semi-global block matching algorithm for non-rectified images
*/

//STL
#include <vector>
#include <array>
#include <algorithm>
#include <ctime>

//Eigen
#include <Eigen/Eigen>

//OpenCV
#include <opencv2/opencv.hpp>

#include "geometry.h"
#include "eucm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"

using namespace std;

using cv::Rect;
using cv::Point;
using cv::Mat_;
using cv::Size;

void EnhancedStereo::computeEpipole()
{
    Vector3d t21 = Transform12.transInv();
    cam2.projectPoint(t21, epipole);
    epipolePx.x = round(epipole[0]);
    epipolePx.y = round(epipole[1]);
}

void EnhancedStereo::computeReconstructed()
{
    vector<Eigen::Vector2d> imagePointVec;
    imagePointVec.reserve(cam1.width*cam1.height);
    for (int v = 0; v < cam1.height; v++)
    {
        for (int u = 0; u < cam1.width; u++)
        {
            imagePointVec.emplace_back(u, v);
        }
    }
    cam1.reconstructPointCloud(imagePointVec, reconstVec);
}

void EnhancedStereo::computeRotated()
{
    Transform12.inverseRotate(reconstVec, reconstRotVec);
}

void EnhancedStereo::computePinf()
{
    cam2.projectPointCloud(reconstRotVec, pinfVec);
    pinfPxVec.resize(pinfVec.size());
    for (int i = 0; i < pinfVec.size(); i++)
    {
        pinfPxVec[i].x = round(pinfVec[i][0]);
        pinfPxVec[i].y = round(pinfVec[i][1]);
    }
        
}

void EnhancedStereo::computeEpipolarCurves()
{
    Vector3d t21 = Transform12.transInv();
    epipolarVec.clear();
    const double & alpha = cam2.params[0];
    const double & beta = cam2.params[1];
    const double & fu = cam2.params[2];
    const double & fv = cam2.params[3];
    const double & u0 = cam2.params[4];
    const double & v0 = cam2.params[5];
    
    double gamma = 1 - alpha;
    double ag = (alpha - gamma);
    double a2b = alpha*alpha*beta;
    
    double fufv = fu * fv;
    double fufu = fu * fu;
    double fvfv = fv * fv;
    
    for (const auto & pt : reconstRotVec)
    {
        epipolarVec.emplace_back();
        Polynomial2 & surf = epipolarVec.back();
        Vector3d plane = pt.cross(t21);
        const double & A = plane[0];
        const double & B = plane[1];
        const double & C = plane[2];
        double AA = A * A;
        double BB = B * B;
        double CC = C * C;
        double CCfufv = CC * fufv;
        if (CCfufv/(AA + BB) < 0.5) // the curve passes through the projection center
        {
            surf.kuu = surf.kuv = surf.kvv = 0;
            surf.ku = A/fu;
            surf.kv = B/fv;
            surf.k1 = -u0*A/fu - v0*B/fv;
        }
        else
        {
            // compute first 4 coefficients directly
            surf.kuu = (AA*ag + CC*a2b)/(CC*fufu);  // kuu
            surf.kuv = 2*A*B*ag/(CCfufv);  // kuv
            surf.kvv = (BB*ag + CC*a2b)/(CC*fvfv);  // kvv
            surf.ku = 2*(-(AA*fv*u0 + A*B*fu*v0)*ag - 
                            A*C*fufv*gamma - CC*a2b*fv*u0)/(CCfufv*fu);  // kv
            surf.kv = 2*(-(BB*fu*v0 + A*B*fv*u0)*ag - 
                            B*C*fufv*gamma - CC*a2b*fu*v0)/(CCfufv*fv);  // kv
                            
            // the last one is computed using the fact that
            // the epipolar curves pass through the epipole
            surf.k1 = -(surf.kuu*epipole[0]*epipole[0] 
                            + surf.kuv*epipole[0]*epipole[1] 
                            + surf.kvv*epipole[1]*epipole[1] 
                            + surf.ku*epipole[0] + surf.kv*epipole[1]);
        }
    }    
}

void EnhancedStereo::traceEpipolarLine(Point pt, Mat_<uint8_t> & out)
{
    int idx = pt.y * cam1.width + pt.x;
    cout << idx << endl;
    cout << pinfPxVec.size() << endl;
    Point pinf = pinfPxVec[idx];
    cout << pinf << endl;
    CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
    
    int count = abs(pinf.x - epipolePx.x);
    for (int i = 0; i < count; i++)
    {
        out(raster.y, raster.x) = 0;
        out(raster.y+1, raster.x) = 0;
        out(raster.y, raster.x+1) = 0;
        out(raster.y+1, raster.x+1) = 0;
        raster.makeStep();
    }
}

void EnhancedStereo::createBuffer()
{
    cout << "create" << endl;
    int bufferWidth = smallWidth()*dispMax;
    if (errorBuffer.cols != bufferWidth or errorBuffer.rows != smallHeight())
    {
        errorBuffer = Mat_<uint8_t>(Size(bufferWidth, smallHeight()));
    }
    if (tableauLeft.cols != bufferWidth or tableauLeft.rows != smallHeight())
    {
        tableauLeft = Mat_<int>(Size(bufferWidth, smallHeight()));
    }
    if (tableauRight.cols != bufferWidth or tableauRight.rows != smallHeight())
    {
        tableauRight = Mat_<int>(Size(bufferWidth, smallHeight()));
    }
    if (tableauTop.cols != bufferWidth or tableauTop.rows != smallHeight())
    {
        tableauTop = Mat_<int>(Size(bufferWidth, smallHeight()));
    }
    if (tableauBottom.cols != bufferWidth or tableauBottom.rows != smallHeight())
    {
        tableauBottom = Mat_<int>(Size(bufferWidth, smallHeight()));
    }
    if (smallDisparity.cols != smallWidth() or smallDisparity.rows != smallHeight())
    {
        smallDisparity = Mat_<uint8_t>(Size(smallWidth(), smallHeight()));
    }
}

void EnhancedStereo::comuteStereo(const Mat_<uint8_t> & img1, 
        const Mat_<uint8_t> & img2,
        Mat_<uint8_t> & disparity)
{
    computeCost(img1, img2);
    computeDynamicProgramming();
    reconstructDisparity();
    upsampleDisparity(img1, disparity);
}

void EnhancedStereo::computeCost(const Mat_<uint8_t> & img1, const Mat_<uint8_t> & img2)
{
    double T1 = 0, T2 = 0; // time profiling
    int halfBlockSize = blockSize / 2;
    cout << "cost" << endl;
    for (int v = 0; v < smallHeight(); v++)
    {
        for (int u = 0; u < smallWidth(); u++)
        {
            int idx = getLinearIdx(vBig(v), uBig(u));
            uint8_t * outPtr = errorBuffer.row(v).data + u*dispMax;
            
            Point pinf = pinfPxVec[idx];
            CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
            
            // the remap Mat
            Mat_<uint8_t> img2remap(Size(blockSize - 1 + dispMax, blockSize));
            
            for (int i = 0; i < img2remap.cols; i++, raster.makeStep())
            {
                for (int j = -halfBlockSize; j <= halfBlockSize; j++)
                {
                    img2remap(halfBlockSize + j, i) = img2(raster.y + j, raster.x + 2);
                }
            }
            
            // compute the actual error
            for (int i = 0; i < dispMax; i++, outPtr++)
            {
                int acc = 0;
                for (int x2 = -halfBlockSize; x2 <= halfBlockSize; x2++)
                {
                    for (int x1 = -halfBlockSize; x1 <= halfBlockSize; x1++)
                    {
                        acc += abs(img1(vBig(v) + x2, uBig(u)- x1) - 
                            img2remap(halfBlockSize + x2, i + halfBlockSize + x1));
                    }
                }
                
                *outPtr = acc / blockSize / blockSize;
            }
        }
    }
    cout << "read " << T1 / CLOCKS_PER_SEC << endl;
    cout << "write " << T2 / CLOCKS_PER_SEC << endl;
} 

void EnhancedStereo::computeDynamicStep(const int* inCost, const uint8_t * error, int * outCost)
{
    int bestCost = inCost[0];
    for (int i = 1; i < dispMax; i++)
    {
        bestCost = min(bestCost, inCost[i]);
    }
    int & val0 = outCost[0];
    val0 = inCost[0];
    val0 = min(val0, inCost[1] + lambdaStep);
    val0 = min(val0, bestCost + lambdaJump);
    val0 += error[0];
    for (int i = 1; i < dispMax-1; i++)
    {
        int & val = outCost[i];
        val = inCost[i];
        val = min(val, inCost[i + 1] + lambdaStep);
        val = min(val, inCost[i - 1] + lambdaStep);
        val = min(val, bestCost + lambdaJump);
        val += error[i];
    }
    int & vald = outCost[dispMax - 1];
    vald = inCost[dispMax - 1];
    vald = min(vald, inCost[dispMax - 2] + lambdaStep);
    vald = min(vald, bestCost + lambdaJump);
    vald += error[dispMax - 1];
}

void EnhancedStereo::computeDynamicProgramming()
{
    cout << "left" << endl;
    // left tableau init
    cout << tableauLeft.size() << endl;
    for (int v = 0; v < smallHeight(); v++)
    {
        int * tableauRow = (int *)(tableauLeft.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        // init the first row
        copy(errorRow, errorRow + dispMax, tableauRow);
        // fill up the tableau
        for (int u = 1; u < smallWidth(); u++)
        {
            computeDynamicStep(tableauRow + (u - 1)*dispMax, errorRow + u*dispMax, tableauRow + u*dispMax);
        }
    }
    cout << "right" << endl;    
    // right tableau init
    for (int v = 0; v < smallHeight(); v++)
    {
        int * tableauRow = (int *)(tableauLeft.row(v).data);
        uint8_t * errorRow = errorBuffer.row(v).data;
        int base = (smallWidth() - 1) * dispMax;
        copy(errorRow + base, errorRow + base + dispMax, tableauRow + base);
        
        for (int u = smallWidth() - 2; u >= 0; u--)
        {
            computeDynamicStep(tableauRow + (u + 1)*dispMax, errorRow + u*dispMax, tableauRow + u*dispMax);
        }
    }
    cout << "top" << endl;
    // top-down tableau init
    for (int u = 0; u < smallWidth(); u++)
    {
        auto tableauCol = tableauTop(Rect(u*dispMax, 0, dispMax, smallHeight()));
        auto errorCol = errorBuffer(Rect(u*dispMax, 0, dispMax, smallHeight()));
        copy(errorCol.data, errorCol.data + dispMax, (int*)(tableauCol.data));
        for (int v = 1; v < smallHeight(); v++)
        {
            computeDynamicStep((int*)(tableauCol.row(v-1).data), 
                    errorCol.row(v).data,
                    (int*)(tableauCol.row(v).data));
        }
    }
    cout << "bottom" << endl;
    // bottom-up tableau init
    for (int u = 0; u < smallWidth(); u++)
    {
        auto tableauCol = tableauBottom(Rect(u*dispMax, 0, dispMax, smallHeight()));
        auto errorCol = errorBuffer(Rect(u*dispMax, 0, dispMax, smallHeight()));
        int vLast = smallHeight() - 1;
        copy(errorCol.row(vLast).data, 
                errorCol.row(vLast).data + dispMax, 
                (int*)(tableauCol.row(vLast).data));
        for (int v = smallHeight() - 2; v >= 0; v--)
        {
            computeDynamicStep((int*)(tableauCol.row(v+1).data), 
                    errorCol.row(v).data,
                    (int*)(tableauCol.row(v).data));
        }
    }
    
}

void EnhancedStereo::reconstructDisparity()
{
    
    for (int v = 0; v < smallHeight(); v++)
    {
        int* dynRow1 = (int*)(tableauLeft.row(v).data);
        int* dynRow2 = (int*)(tableauRight.row(v).data);
        int* dynRow3 = (int*)(tableauTop.row(v).data);
        int* dynRow4 = (int*)(tableauBottom.row(v).data);
        for (int u = 0; u < smallWidth(); u++)
        {
            int bestCost = 10000;
            uint8_t & bestDisp = smallDisparity(v, u);
            bestDisp = 0;
            for (int d = 0; d < dispMax; d++)
            {
                int base = u * dispMax;
                int acc = dynRow1[base + d] + dynRow2[base + d] + dynRow3[base + d] + dynRow4[base + d];
                if (bestCost > acc)
                {
                    bestCost = acc;
                    bestDisp = d;
                }
            }
        }
    }
    medianBlur(smallDisparity, smallDisparity, 3);
}


void EnhancedStereo::upsampleDisparity(const Mat_<uint8_t> & img1, Mat_<uint8_t> & disparity)
{
    resize(smallDisparity, disparity, Size(0, 0), blockSize, blockSize, 0);
}

