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
Stereo vision definition.
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

void StereoEUCM::computeEpipole()
{
    Vector3d t21 = Transform12.transInv();
    cam2.projectPoint(t21, epipole);
    epipolePx.x = round(epipole[0]);
    epipolePx.y = round(epipole[1]);
}

void StereoEUCM::computeReconstructed()
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

void StereoEUCM::computeRotated()
{
    Transform12.inverseRotate(reconstVec, reconstRotVec);
}

void StereoEUCM::computePinf()
{
    cam2.projectPointCloud(reconstRotVec, pinfVec);
    pinfPxVec.resize(pinfVec.size());
    for (int i = 0; i < pinfVec.size(); i++)
    {
        pinfPxVec[i].x = round(pinfVec[i][0]);
        pinfPxVec[i].y = round(pinfVec[i][1]);
    }
        
}

void StereoEUCM::computeEpipolarCurves()
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

void StereoEUCM::traceEpipolarLine(Point pt, Mat_<uint8_t> & out)
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

void StereoEUCM::computeCost(const cv::Mat_<uint8_t> & img1,
        const cv::Mat_<uint8_t> & img2, cv::Mat_<uint8_t> & out)
{
    assert(img1.cols == cam1.width and img1.rows == cam1.height);
    int bufferWidth = floor(img1.cols / blockSize)*dispMax;
    int bufferHeight = floor(img1.rows / blockSize);
    if (bufferWidth != out.cols or bufferHeight != out.rows)
    {
        out = Mat_<uint8_t>(Size(bufferWidth, bufferHeight));
    }
    double T1 = 0, T2 = 0;
    int halfBlockSize = blockSize / 2;
    for (int v = verticalMargin; v < img1.rows - verticalMargin; v += blockSize)
    {
        for (int u = 2*dispMax; u < img1.cols - 2*dispMax; u += blockSize)
        {
            int idx = cam1.width * v + u;
            uint8_t * outPtr = out.data + int(floor(u / blockSize)*dispMax + bufferWidth*floor(v / blockSize));
            
            Point pinf = pinfPxVec[idx];
            CurveRasterizer<Polynomial2> raster(pinf.x, pinf.y, epipolePx.x, epipolePx.y, epipolarVec[idx]);
            
            // the remap Mat
            Mat_<uint8_t> img2remap(Size(blockSize - 1 + dispMax, blockSize));
            
//            auto t1 = clock();
            for (int i = 0; i < img2remap.cols; i++, raster.makeStep())
            {
                for (int j = -halfBlockSize; j <= halfBlockSize; j++)
                {
                    img2remap(halfBlockSize + j, i) = img2(raster.y + j, raster.x + 2);
                }
            }
            
//            cv::imshow("remap", img2remap);
//            cv::waitKey();
//            auto t2 = clock();
//            T1 += double(t2 - t1);
//            t1 = clock();
            for (int i = 0; i < dispMax; i++, outPtr++)
            {
                int acc = 0;
                for (int x2 = -halfBlockSize; x2 <= halfBlockSize; x2++)
                {
                    for (int x1 = -halfBlockSize; x1 <= halfBlockSize; x1++)
                    {
                        acc += abs(img1(v + x2, u - x1) - 
                            img2remap(halfBlockSize + x2, i + halfBlockSize + x1));
                    }
                }
                
                *outPtr = acc / blockSize / blockSize;
//                cout << int(*outPtr) << endl;
            }
//            cout << endl;
//            t2 = clock();
//            T2 += double(t2 - t1);
        }
    }
    cout << "read " << T1 / CLOCKS_PER_SEC << endl;
    cout << "write " << T2 / CLOCKS_PER_SEC << endl;
} 


void StereoEUCM::computeDynamicProgramming(const cv::Mat_<uint8_t> & costMat, cv::Mat_<int> & out)
{
    if (out.cols != costMat.cols or out.rows != costMat.rows*3)
    {
        // out contains 4 tableaux for DP
        // the first two are horizontal (from left, from right)
        // the second two are vertical (top-down, bottom-up);
        out = Mat_<int>(costMat.rows*4, costMat.cols);
        cout << out.size() << endl;
        out.setTo(0);
        cout << "Reset" << endl;
    }
    cout << "out created " << out.size() << endl;
    cout << costMat.rows << endl;
    int lambdaJump = 32;
    int lambdaStep = 8;
    
    // left tableau init
    for (int v = verticalMargin/blockSize; v < costMat.rows - verticalMargin/blockSize; v++)
    {
        int bestCost = 256;
        auto outRow = out.row(v);
        auto costRow = costMat.row(v);
        for (int d = 0; d < dispMax; d++)
        {
            outRow(d) = costRow(d);
            bestCost = min( bestCost, int(outRow(d)) );
        }
        for (int u = 1; u < costMat.cols / dispMax; u++)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outRow((u - 1)*dispMax);
            val0 = min(val0, outRow((u - 1)*dispMax + 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outRow(u*dispMax) = val0 + costRow(u*dispMax);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outRow((u - 1)*dispMax + d);
                val = min(val, outRow((u - 1)*dispMax + d + 1) + lambdaStep);
                val = min(val, outRow((u - 1)*dispMax + d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outRow(u*dispMax + d) = val + costRow(u*dispMax + d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
            int vald = outRow(u*dispMax - 1);
            vald = min(vald, outRow(u*dispMax - 2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outRow((u + 1)*dispMax - 1) = vald + costRow((u + 1)*dispMax - 1);
            bestCost = outRow(u*dispMax);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outRow(u*dispMax + d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
        
    // right tableau init
    for (int v =  verticalMargin/blockSize; v < costMat.rows - verticalMargin/blockSize; v++)
    {
//        cout << v << endl;
        int bestCost = 256;
        auto outRow = out.row(costMat.rows + v);
        auto costRow = costMat.row(v);
        for (int u = costMat.cols - dispMax; u < costMat.cols; u++)
        {
            outRow(u) = costRow(u);
            bestCost = min( bestCost, int(outRow(u)) );
        }
//        cout << 111 << endl;
        for (int u = costMat.cols / dispMax - 2; u >= 0; u--)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outRow((u + 1)*dispMax);
            val0 = min(val0, outRow((u + 1)*dispMax + 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outRow(u*dispMax) = val0 + costRow(u*dispMax);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outRow((u + 1)*dispMax + d);
                val = min(val, outRow((u + 1)*dispMax  + d + 1) + lambdaStep);
                val = min(val, outRow((u + 1)*dispMax  + d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outRow(u*dispMax + d) = val + costRow(u*dispMax + d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outRow((u + 2)*dispMax - 1);
            vald = min(vald, outRow((u + 2)*dispMax - 2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outRow((u + 1)*dispMax - 1) = vald + costRow((u + 1)*dispMax - 1);
            bestCost = outRow(u*dispMax);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outRow(u*dispMax + d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
//    lambdaJump = 100;
//    lambdaStep = 1;
    const int G = 2;
    // top-down tableau init
    for (int x =  0; x < costMat.cols; x+= dispMax)
    {
        int bestCost = 256;
        auto outMat = out(Rect(x, 2*costMat.rows, dispMax, costMat.rows));
        auto costSeg1 = out(Rect(x, 0, dispMax, costMat.rows));
        auto costSeg2 = out(Rect(x, costMat.rows, dispMax, costMat.rows));
        for (int u = 0; u < dispMax; u++)
        {
            outMat(0, u) = (costSeg1(0, u) + costSeg2(0, u))/G;
            bestCost = min( bestCost, int(outMat(0, u)) );
        }
//        cout << 111 << endl;
        for (int v = 1; v < costMat.rows; v++)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outMat(v - 1, 0);
            val0 = min(val0, outMat(v - 1, 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outMat(v, 0) = val0 + (costSeg1(v, 0) + costSeg2(v, 0))/G;
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outMat(v - 1, d);
                val = min(val, outMat(v - 1, d + 1) + lambdaStep);
                val = min(val, outMat(v - 1, d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outMat(v, d) = val + (costSeg1(v, d) + costSeg2(v, d))/G;
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outMat(v - 1, dispMax-1);
            vald = min(vald, outMat(v - 1, dispMax-2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outMat(v, dispMax-1) = vald + (costSeg1(v, dispMax-1) + costSeg2(v, dispMax-1))/G;
            bestCost = outMat(v, 0);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outMat(v, d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
    // bottom-up tableau init
    for (int x =  0; x < costMat.cols; x+= dispMax)
    {
        int bestCost = 256;
        auto outMat = out(Rect(x, 2*costMat.rows, dispMax, costMat.rows));
        auto costSeg1 = out(Rect(x, 0, dispMax, costMat.rows));
        auto costSeg2 = out(Rect(x, costMat.rows, dispMax, costMat.rows));
        for (int u = 0; u < dispMax; u++)
        {
            outMat(costMat.rows-1, u) = (costSeg1(costMat.rows-1, u) + costSeg2(costMat.rows-1, u))/G;
            bestCost = min( bestCost, int(outMat(costMat.rows, u)) );
        }
//        cout << 111 << endl;
        for (int v = costMat.rows - 2; v >= 0; v--)
        {
//            cout << "COL : " << v << endl;
            int val0;
            val0 = outMat(v + 1, 0);
            val0 = min(val0, outMat(v + 1, 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outMat(v, 0) = val0 + (costSeg1(v, 0) + costSeg2(v, 0))/G;
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outMat(v + 1, d);
                val = min(val, outMat(v + 1, d + 1) + lambdaStep);
                val = min(val, outMat(v + 1, d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outMat(v, d) = val + (costSeg1(v, d)  + costSeg2(v, d))/G;
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outMat(v + 1, dispMax-1);
            vald = min(vald, outMat(v + 1, dispMax-2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outMat(v, dispMax-1) = vald + (costSeg1(v, dispMax-1) + costSeg2(v, dispMax-1))/G;
            bestCost = outMat(v, 0);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outMat(v, d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
    //left tableau computation
}


/*
void StereoEUCM::computeDynamicProgramming(const cv::Mat_<uint8_t> & costMat, cv::Mat_<int> & out)
{
    if (out.cols != costMat.cols or out.rows != costMat.rows*3)
    {
        // out contains 4 tableaux for DP
        // the first two are horizontal (from left, from right)
        // the second two are vertical (top-down, bottom-up);
        out = Mat_<int>(costMat.rows*4, costMat.cols);
        cout << out.size() << endl;
        out.setTo(0);
        cout << "Reset" << endl;
    }
    cout << "out created " << out.size() << endl;
    cout << costMat.rows << endl;
    int lambdaJump = 32;
    int lambdaStep = 8;
    
    // left tableau init
    for (int v = verticalMargin/blockSize; v < costMat.rows - verticalMargin/blockSize; v++)
    {
        int bestCost = 256;
        auto outRow = out.row(v);
        auto costRow = costMat.row(v);
        for (int d = 0; d < dispMax; d++)
        {
            outRow(d) = costRow(d);
            bestCost = min( bestCost, int(outRow(d)) );
        }
        for (int u = 1; u < costMat.cols / dispMax; u++)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outRow((u - 1)*dispMax);
            val0 = min(val0, outRow((u - 1)*dispMax + 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outRow(u*dispMax) = val0 + costRow(u*dispMax);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outRow((u - 1)*dispMax + d);
                val = min(val, outRow((u - 1)*dispMax + d + 1) + lambdaStep);
                val = min(val, outRow((u - 1)*dispMax + d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outRow(u*dispMax + d) = val + costRow(u*dispMax + d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
            int vald = outRow(u*dispMax - 1);
            vald = min(vald, outRow(u*dispMax - 2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outRow((u + 1)*dispMax - 1) = vald + costRow((u + 1)*dispMax - 1);
            bestCost = outRow(u*dispMax);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outRow(u*dispMax + d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
        
    // right tableau init
    for (int v =  verticalMargin/blockSize; v < costMat.rows - verticalMargin/blockSize; v++)
    {
//        cout << v << endl;
        int bestCost = 256;
        auto outRow = out.row(costMat.rows + v);
        auto costRow = costMat.row(v);
        for (int u = costMat.cols - dispMax; u < costMat.cols; u++)
        {
            outRow(u) = costRow(u);
            bestCost = min( bestCost, int(outRow(u)) );
        }
//        cout << 111 << endl;
        for (int u = costMat.cols / dispMax - 2; u >= 0; u--)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outRow((u + 1)*dispMax);
            val0 = min(val0, outRow((u + 1)*dispMax + 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outRow(u*dispMax) = val0 + costRow(u*dispMax);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outRow((u + 1)*dispMax + d);
                val = min(val, outRow((u + 1)*dispMax  + d + 1) + lambdaStep);
                val = min(val, outRow((u + 1)*dispMax  + d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outRow(u*dispMax + d) = val + costRow(u*dispMax + d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outRow((u + 2)*dispMax - 1);
            vald = min(vald, outRow((u + 2)*dispMax - 2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outRow((u + 1)*dispMax - 1) = vald + costRow((u + 1)*dispMax - 1);
            bestCost = outRow(u*dispMax);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outRow(u*dispMax + d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
    // top-down tableau init
    for (int x =  0; x < costMat.cols; x+= dispMax)
    {
        int bestCost = 256;
        auto outMat = out(Rect(x, 2*costMat.rows, dispMax, costMat.rows));
        auto costSeg = costMat(Rect(x, 0, dispMax, costMat.rows));
        for (int u = 0; u < dispMax; u++)
        {
            outMat(0, u) = costSeg(0, u);
            bestCost = min( bestCost, int(outMat(0, u)) );
        }
//        cout << 111 << endl;
        for (int v = 1; v < costMat.rows; v++)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outMat(v - 1, 0);
            val0 = min(val0, outMat(v - 1, 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outMat(v, 0) = val0 + costSeg(v, 0);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outMat(v - 1, d);
                val = min(val, outMat(v - 1, d + 1) + lambdaStep);
                val = min(val, outMat(v - 1, d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outMat(v, d) = val + costSeg(v, d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outMat(v - 1, dispMax-1);
            vald = min(vald, outMat(v - 1, dispMax-2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outMat(v, dispMax-1) = vald + costSeg(v, dispMax-1);
            bestCost = outMat(v, 0);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outMat(v, d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
    // bottom-up tableau init
    for (int x =  0; x < costMat.cols; x+= dispMax)
    {
        int bestCost = 256;
        auto outMat = out(Rect(x, 2*costMat.rows, dispMax, costMat.rows));
        auto costSeg = costMat(Rect(x, 0, dispMax, costMat.rows));
        for (int u = 0; u < dispMax; u++)
        {
            outMat(costMat.rows-1, u) = costSeg(costMat.rows-1, u);
            bestCost = min( bestCost, int(outMat(costMat.rows, u)) );
        }
//        cout << 111 << endl;
        for (int v = costMat.rows - 2; v >= 0; v--)
        {
//            cout << "COL : " << u << endl;
            int val0;
            val0 = outMat(v + 1, 0);
            val0 = min(val0, outMat(v + 1, 1) + lambdaStep);
            val0 = min(val0, bestCost + lambdaJump);
            outMat(v, 0) = val0 + costSeg(v, 0);
            for (int d = 1; d < dispMax-1; d++)
            {
                int val = outMat(v + 1, d);
                val = min(val, outMat(v + 1, d + 1) + lambdaStep);
                val = min(val, outMat(v + 1, d - 1) + lambdaStep);
                val = min(val, bestCost + lambdaJump);
                outMat(v, d) = val + costSeg(v, d);
//                cout << d << " " << val + costRow(u*dispMax + d) << " " 
//              << int(costRow(u*dispMax + d)) << endl;
            }
//            cout << u << endl;
            int vald = outMat(v + 1, dispMax-1);
            vald = min(vald, outMat(v + 1, dispMax-2) + lambdaStep);
            vald = min(vald, bestCost + lambdaJump);
            outMat(v, dispMax-1) = vald + costSeg(v, dispMax-1);
            bestCost = outMat(v, 0);
            for (int d = 1; d < dispMax; d++)
            {
                bestCost = min(bestCost, int(outMat(v, d)));
            }
//            cout << "Best : " << bestCost << endl;
//            cout << endl;
        }
        
    }
    
    //left tableau computation
}*/

/*
class StereoEUCM
{
public:
    StereoEUCM(Transformation<double> T12, int imageWidth, int imageHeight,
            const array<double, 6> & params1, const array<double, 6> & params2)
            : Transform12(T12), 
            cam1(imageWidth, imageHeight, params1.data()), 
            cam2(imageWidth, imageHeight, params2.data()) {}

    ~StereoEUCM();

    void setTransformation(Transformation<double> T12) { Transform12 = T12; } 
    
    // f2(t21) -- projection of the first camera's projection center
    void computeEpipole();
    
    // computes reconstVec -- reconstruction of every pixel of the first image
    void computeReconstructed();
    
    // computes reconstRotVec -- reconstVec rotated into the second frame
    void computeRotated();
    
    // computes pinfVec -- projections of all the reconstructed points from the first image
    // onto the second image as if they were at infinity
    void computePinf();
    
    // calculate the coefficients of the polynomials for all the 
    void computeEpipolarCurves();
    
    // fills up the output with photometric errors between the val = I1(pi) and 
    // the values from I2 on the epipolar curve
    void computeError(EpipolarRasterizer rasterizer, uint8_t val, cv::Mat_<uint8_t> img2, uint8_t * output);
    
private:
    Transformation<double> Transform12;  // pose of the first to the second camera
    EUCM cam1, cam2;
   
    vector<Eigen::Vector3d> reconstVec;  // reconstruction of every pixel by cam1
    vector<Eigen::Vector3d> reconstRotVec;  // reconstVec rotated into the second frame
    
    Eigen::Vector2d epipole;
    vector<Eigen::Vector2d> pinfVec;  // projection of reconstRotVec by cam2
    
    cv::Point epipolePx;  // Px means that it is discrete
    vector<cv::Point> pinfPxVec;
    
    vector<Polynomial2> epipolarVec;
};
*/


