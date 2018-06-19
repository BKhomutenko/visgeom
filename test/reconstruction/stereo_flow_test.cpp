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

#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "json.h"

#include "geometry/geometry.h"
#include "projection/eucm.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "render/render.h"

#include <opencv2/highgui/highgui.hpp>  // Video write
//FIXME make an argument
ofstream results;
string histDataName;
void analyzeError(const Mat32f & depthGT, Mat32f & depth, 
        const Mat32f & sigma, const ScaleParameters & scaleParams)
{
    Mat8u inlierMat(depth.size());
    inlierMat.setTo(0);
    int Nmax = 0, Ngt = 0;
    double dist = 0;
    int N = 0;
    double err = 0, err2 = 0;
    double err3 = 0;
    std::ofstream ofs;
  ofs.open (histDataName, std::ofstream::out | std::ofstream::app);


 
    for (int u = 0; u < depth.cols; u++)
    {
        for (int v = 0; v < depth.rows; v++)
        {
            int ugt = scaleParams.uConv(u);
            int vgt = scaleParams.vConv(v);
            if (depthGT(vgt, ugt) != 0) Ngt++;
            
            if (depthGT(vgt, ugt) == 0 or depth(v, u) == 0 )
            {
                depth(v, u) = 0;
                continue;
                
            }
            
            if (depthGT(vgt, ugt) != depthGT(vgt, ugt) or depth(v, u) != depth(v, u)) continue;
            Nmax++;
            dist += depthGT(vgt, ugt);
            ofs << (depthGT(vgt, ugt) - depth(v, u)) << "    " << depthGT(vgt, ugt) << "   " << sigma(v, u) << endl;
            if (sigma(v, u) > 10 or abs(depthGT(vgt, ugt) - depth(v, u)) >  sigma(v, u))
            {
                continue;
            }
            inlierMat(v, u) = 255;
            err += depthGT(vgt, ugt) - depth(v, u);
            err2 += pow(depthGT(vgt, ugt) - depth(v, u), 2);
            N++;
        }
    }
    cout << "avg err : " << err / N *1000 << " avg err2 : " 
        << sqrt(err2 / N)*1000  << " number of inliers : " << 100 * N / double(Nmax)
        << "   average distance : " << dist / Nmax << endl;
    
    results << sqrt(err2 / N)*1000  << "    " << 100 * N / double(Nmax)
        << "    " << dist / Nmax << "    " << N / double(Ngt);
     ofs.close();    
    imshow("inliers", inlierMat);
}

int main(int argc, char** argv) 
{

    ptree root;
    read_json(argv[1], root);
    
    Transf xiCam12 = readTransform(root.get_child("Extrinsics"));
    
    EnhancedCamera camera1( readVector<double>(root.get_child("intrinsics1")).data() );
    EnhancedCamera camera2( readVector<double>(root.get_child("intrinsics2")).data() );
    
    //init stereoParameters
    SgmParameters stereoParams(root.get_child("stereo_parameters"));
    
    //the stereo is computed backwards
    EnhancedSgm sgm(xiCam12, &camera1, &camera2, stereoParams);
    

    //depth GT
    Mat32f depthGT, depth, sigmaMat;
    Mat8u depth8u;
    DepthMap depthStereo;
    
    //Video writing
    Size S = Size(stereoParams.xMax, 2*stereoParams.yMax);

    cv::VideoWriter outputVideo;                                        // Open the output
    outputVideo.open("video.avi", CV_FOURCC('M','J','P','G'), 20, S, true);

    if (!outputVideo.isOpened())
    {
        cout  << "Could not open the output video for write" << endl;
        return -1;
    }
    
    for (auto it1 = root.get_child("img1").begin(), it2 = root.get_child("img2").begin();
         it1 != root.get_child("img1").end(); ++it1, ++it2)
    {

        Mat8u img1 = imread(it1->second.get_value<string>(), 0);
        Mat8u img2 = imread(it2->second.get_value<string>(), 0);
        sgm.computeStereo(img1, img2, depthStereo);
       
        
        depthStereo.toInverseMat(depth);
        depth *= 800;
        depth.copyTo(depth8u);
//        depthStereo.toMat(depth);
//        resize(depth8u, depth8u, Size(0,0), 3, 3, cv::INTER_NEAREST);
        
        Mat8uc3 heatMap;
        cv::cvtColor(depth8u, heatMap, CV_GRAY2BGR);
        cv::applyColorMap(heatMap, heatMap, cv::COLORMAP_JET);
        
        for (int v = 0; v < depth8u.rows; v++)
        {
            for (int u = 0; u < depth8u.cols; u++)
            {
                if (depth8u(v, u) == 0) heatMap(v, u) = cv::Vec<uchar, 3>(0, 0, 0);
                
            }
        }
        
        imshow("img1", img1);
        imshow("img2", img2);
        imshow("depth", heatMap);
        
        Mat imgCut = img1.rowRange(stereoParams.v0,
                stereoParams.v0 + stereoParams.yMax).colRange(stereoParams.u0, 
                stereoParams.u0 + stereoParams.xMax);
        Mat8uc3 res(heatMap.rows*2, heatMap.cols);
        Mat8uc3 imgCutColor;
        cv::cvtColor(imgCut, imgCutColor, CV_GRAY2BGR);
        imgCutColor.copyTo(res.rowRange(0, imgCutColor.rows));
        heatMap.copyTo(res.rowRange(imgCutColor.rows, 2*imgCutColor.rows));
        cout << res.rowRange(0, imgCutColor.rows).size() << endl;
        cout << imgCutColor.size() << endl;
        outputVideo << res;
        waitKey(30);
    }
    results.close();
    return 0;
}



