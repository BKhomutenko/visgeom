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


#include "reconstruction/eucm_motion_stereo.h"

#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "utils/filter.h"
#include "geometry/geometry.h"
#include "projection/eucm.h"

#include "reconstruction/stereo_misc.h"
#include "reconstruction/eucm_epipolar.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/epipolar_descriptor.h"

//FIXME the same as DepthMap::filter
void filter(double & v1, double & s1, const double v2, const double s2)
{
    double K = 1. / (s1 + s2);
    v1 = (v1 * s2 + v2 * s1) * K;
    s1 = max(s1 * s2 * K, 0.1);
}

bool MotionStereo::selectPoint(int x, int y)
{
    gu = _params.uConv(x);
    gv = _params.vConv(y);
    flags |= GLB_UV;
    
    //--Check point's saliency
    if (_maskMat(gv, gu) < _params.gradientThresh) return false;
    
    Vector2d pt(gu, gv);
    
    if (not _camera1->reconstructPoint(pt, gX)) return false;
    flags |= GLB_X;
    
    Vector2i pti = round(pt);
    auto useInverted = epipoles().chooseEpipole(CAMERA_1, pti, _params.epipoleMargin);
    if (useInverted & EPIPOLE_TOO_CLOSE) return false;
    Vector2i goal = epipoles().getPx(CAMERA_1, useInverted);
    CurveRasterizer<int, Polynomial2> descRaster(round(pt), goal,
                                _epipolarCurves.get(CAMERA_1, gX));
    if (useInverted) descRaster.setStep(-1);

    //to compute one step for the uncertainty estimation
    CurveRasterizer<int, Polynomial2> descRasterUncert = descRaster;
    
    gstep = _epipolarDescriptor.compute(_img1, descRaster, gdescriptor);
    
    descRasterUncert.setStep(gstep);
    descRasterUncert.step();
    
    gu2 = descRasterUncert.u;
    gv2 = descRasterUncert.v;
    
    if (gstep < 1 or not _epipolarDescriptor.goodResp()) return false;
    flags |= GLB_STEP | GLB_DESCRIPTOR;
    return true;
}


bool MotionStereo::computeUncertainty(double d, double s)
{
    //TODO replace assert?
    uint32_t neededFlag = GLB_X;
    assert( (flags & neededFlag) ^ neededFlag == 0);
    
    
    if (d == OUT_OF_RANGE)  // no prior
    {
        //just rotate
//        return false; //FIXME
        Vector3d Xmax = R21() * gX;
        if (not _camera2->projectPoint(Xmax, gptStart)) return false;
        gptStartRound = round(gptStart);
        auto useInverted = epipoles().chooseEpipole(CAMERA_2, gptStartRound, _params.epipoleMargin);
        if (useInverted & EPIPOLE_TOO_CLOSE) return false;
        gptFinRound = epipoles().getPx(CAMERA_2, useInverted);
        if (useInverted & EPIPOLE_INVERTED)
        {
            gdispMax = _params.dispMax;
            flags |= GLB_INVERTED_SAMPLING;
        }
        else
        {
            int delta = round( max( abs(gptStartRound[0] - gptFinRound[0]),
                                     abs(gptStartRound[1] - gptFinRound[1]) ) );
            gdispMax = min(_params.dispMax, delta);
        }
        
        flags |= GLB_START_POINT | GLB_DISP_MAX;
    }
    else // there is a prior
    {
        gX.normalize();
        Vector3d Xmax = gX * (d + 2.5 * s);
        Vector3d Xmin = gX * max(d - 2.5 * s, MIN_DEPTH);
        Xmax = R21() * (Xmax - t12());
        Xmin = R21() * (Xmin - t12());
        Vector2d ptFin;
        if (not _camera2->projectPoint(Xmax, gptStart)) return false;
        if (not _camera2->projectPoint(Xmin, ptFin)) return false;
        int delta = round( max(abs(ptFin[0] - gptStart[0]), abs(ptFin[1] - gptStart[1])) );
        gdispMax = min( _params.dispMax, delta);
        gptStartRound = round(gptStart);
        gptFinRound = round(ptFin);
        flags |= GLB_START_POINT | GLB_DISP_MAX;
    }
    return true;
}

bool MotionStereo::sampleImage(const Mat8u & img2)
{
    uint32_t neededFlag = GLB_START_POINT | GLB_DISP_MAX | GLB_STEP | GLB_X;
    assert(flags & neededFlag == neededFlag);
    
    int distance = gdispMax / gstep + MARGIN;
    
    CurveRasterizer<int, Polynomial2> raster(gptStartRound, gptFinRound,
                                _epipolarCurves.get(CAMERA_2, gX));
    if (flags & GLB_INVERTED_SAMPLING)
    {
        raster.setStep(-1);
    }
    //Important : Epipolar curves are accessed by the reconstructed point in the FIRST frame
                                
    raster.setStep(gstep);
    raster.steps(-HALF_LENGTH);
    
    guVec.clear();
    guVec.reserve(distance);
    gvVec.clear();
    gvVec.reserve(distance);
    gsampleVec.clear();
    gsampleVec.reserve(distance);
    for (int d = 0; d < distance; d++, raster.step())
    {
        if (raster.v < 0 or raster.v >= img2.rows 
            or raster.u < 0 or raster.u >= img2.cols) 
        {
            return false;
        }//sampleVec.push_back(0);
        
        gsampleVec.push_back(img2(raster.v, raster.u));
        guVec.push_back(raster.u);
        gvVec.push_back(raster.v);
    }
    assert(gsampleVec.size() > MARGIN);
    flags |= GLB_SAMPLE_VEC | GLB_UV_VEC;
    return true;
}

void MotionStereo::reconstruct(double & dist, double & sigma, double & cost)
{
    uint32_t neededFlag = GLB_SAMPLE_VEC | GLB_UV | GLB_UV_VEC | GLB_DESCRIPTOR;
    assert(flags & neededFlag == neededFlag);
    
    vector<int> costVec = compareDescriptor(gdescriptor, gsampleVec, _params.flawCost);
    auto bestCostIter = min_element(costVec.begin() + HALF_LENGTH, costVec.end() - HALF_LENGTH);
    
    
    
    if (1/* *bestCostIter < _params.maxError and *bestCostIter < 2*cost FIXME*/)
    {
        int dBest = bestCostIter - costVec.begin();
//        cout << setw(8) << dBest;
        double distNew, sigmaNew;
        triangulate(Vector2d(gu, gv),
                    Vector2d(gu2, gv2),
                    Vector2d(guVec[dBest], gvVec[dBest]), 
                    Vector2d(guVec[dBest + 1], gvVec[dBest + 1]),
                    distNew, sigmaNew); 
        
        if (dist != OUT_OF_RANGE)
        {
            if (abs(distNew - dist) > 2.6*sigma) 
            {
                count_out++;
    //            dist = OUT_OF_RANGE;
    //            cout << gu << " " << gv << " ; " << guVec[dBest] << " " << gvVec[dBest] << " " ;
    //            cout << dist << "+-" << sigma << " ; " << distNew << "+-" << sigmaNew << endl;
            }
    //        else
    //        { count_in++;
            //temporal filtering
            
            filter(dist, sigma, distNew, sigmaNew);
            cost = cost * 0.7 +  (*bestCostIter) * 0.3; //FIXME experimental
        }
        else
        {
//            dist =  dBest;            
            dist = distNew;
            sigma = sigmaNew;
            cost = *bestCostIter;
        }    
////        const double K = 1 / (sigma + sigmaNew);
//        dist = (sigma * distNew + sigmaNew * dist) * K;
//        sigma = min(sigma, sigmaNew);  //FIXME an overestimation?
        
        /*
        if (sigma > 1 and gu < 600 and gu > 400 and gv > 300 and gv < 550 )
        {
            cout << sigma << " " << dist << endl;
            cout    << gu << " " << gv << " " 
                    <<  guVec[dBest] << " "  << gvVec[dBest] << " " 
                     << guVec[dBest + 1] << " "  << gvVec[dBest + 1] << endl;
        }
        */
//        }
        
    }
    
//    if (gu > 220 and gu < 235  and gv > 103 and gv < 113 )
    /*if (gu > 110 and gu < 126  and gv > 255 and gv < 265 )
    {
        cout << gu << "   " << gv << endl;
        cout << "depth: " << dist
            << " +-" << sigma
            << endl;
        cout << "samples:" << endl;
        for (auto & x : gsampleVec)
        {
            cout << setw(6) << int(x);
        }
        cout << endl;
        cout << "coordinates:" << endl;
        for (auto & x : guVec)
        {
            cout << setw(6) << int(x);
        }
        cout << endl;
        for (auto & x : gvVec)
        {
            cout << setw(6) << int(x);
        }
        cout << endl;
        cout << "descriptor:" << endl;
        for (auto & x : gdescriptor)
        {
            cout << setw(6) << int(x);
        }
        cout << endl;
        cout << "cost:" << endl;
        for (auto & x : costVec)
        {
            cout << setw(6) << int(x);
        }
        cout << endl<< endl;
        
        
    }*/
}

DepthMap MotionStereo::compute(Transf T12, const Mat8u & img2)
{
    //init necessary data structures
    setTransformation(T12);
    DepthMap depthOut(_camera1, _params);
    depthOut.setTo(OUT_OF_RANGE, OUT_OF_RANGE, OUT_OF_RANGE);

//    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    //for each point
    for (int y = 0; y < depthOut.yMax; y++)
    {
        for (int x = 0; x < depthOut.xMax; x++)
        {
//            count4++;
            flags = 0;
            if (not selectPoint(x, y)) 
            { 
//                count1++; 
                continue;
            }
            
            if (not computeUncertainty(OUT_OF_RANGE, OUT_OF_RANGE)) 
            { 
//                count2++; 
                continue;
            }
            
            if (not sampleImage(img2)) {
//                count3++;
                continue;  
            }
            
            reconstruct(depthOut.at(x, y), depthOut.sigma(x, y), depthOut.cost(x, y));
        }
    }
//    cout << count1 << endl;
//    cout << count2 << endl;
//    cout << count3 << endl;
//    cout << count4 << endl;
    
    return depthOut;
}

DepthMap MotionStereo::compute(Transf T12, const Mat8u & img2, const DepthMap & depthIn)
{
    //init necessary data structures
    setTransformation(T12);
    assert(ScaleParameters(depthIn) == ScaleParameters(_params));
    DepthMap depthOut = depthIn;
    
    count_in = 0;
    count_out = 0;
    
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    //for each point
    for (int y = 0; y < depthOut.yMax; y++)
    {
        for (int x = 0; x < depthOut.xMax; x++)
        {
            flags = 0;
            if (not selectPoint(x, y)) { count1++; continue;}
            
            if (not computeUncertainty(depthIn.at(x, y), depthIn.sigma(x, y))) { count2++; continue;}
            
            if (gdispMax / gstep < 2) 
            {
                count3++;
                continue;  
            }  // the uncertainty is too small
            
            // TODO if the uncertainty is small, fuse the two measurements 
            // replace the old one otherwise
            
            if (not sampleImage(img2)) { count4++; continue;}
            
            reconstruct(depthOut.at(x, y), depthOut.sigma(x, y), depthOut.cost(x, y));
        }
    }
//    cout << count1 << endl;
//    cout << count2 << endl;
//    cout << count3 << endl;
//    cout << count4 << endl << endl;
//    cout << count_in << endl;
//    cout << count_out << endl;
    return depthOut;
}


