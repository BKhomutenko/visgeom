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
#include "utils/image_generator.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

//FIXME make an argument
ofstream results;

void analyzeError(const Mat32f & depthGT, const Mat32f & depth, 
        const Mat32f & sigma, const ScaleParameters & scaleParams)
{
    Mat8u inlierMat(depth.size());
    inlierMat.setTo(0);
    int Nmax = 0;
    double dist = 0;
    int N = 0;
    double err = 0, err2 = 0;
    for (int u = 0; u < depth.cols; u++)
    {
        for (int v = 0; v < depth.rows; v++)
        {
            int ugt = scaleParams.uConv(u);
            int vgt = scaleParams.vConv(v);
            if (depthGT(vgt, ugt) == 0 or depth(v, u) == 0 ) continue;
            Nmax++;
            dist += depthGT(vgt, ugt);
            if (depthGT(vgt, ugt) != depthGT(vgt, ugt) or depth(v, u) != depth(v, u)) continue;
            if (sigma(v, u) > 0.3 or abs(depthGT(vgt, ugt) - depth(v, u)) > 2.5 * sigma(v, u)) continue;
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
        << "    " << dist / Nmax;
        
    imshow("inliers", inlierMat);
}

int main(int argc, char** argv) 
{
    ptree root;
    read_json(argv[1], root);
    const vector<double> intrinsic = readVector(root.get_child("camera_intrinsics"));
    const int width = root.get<int>("image.width");
    const int height = root.get<int>("image.height");
    Transf xiCam0 = readTransform(root.get_child("camera_transform"));
    
    Mat8u foreImg = imread(root.get<string>("foreground"), 0);
    
    EnhancedCamera camera(width, height, intrinsic.data());
    
    //init stereoParameters
    SGMParameters stereoParams;
    
    stereoParams.verbosity = root.get<int>("stereo.verbosity");
    stereoParams.salientPoints = true;
    stereoParams.u0 = root.get<int>("stereo.u0");
    stereoParams.v0 = root.get<int>("stereo.v0");
    stereoParams.dispMax = root.get<int>("stereo.disparity_max");
    stereoParams.scale = root.get<int>("stereo.scale");
    stereoParams.flawCost = root.get<int>("stereo.flaw_cost");
    stereoParams.scaleVec = readIntVector(root.get_child("stereo.scale_vector"));
    stereoParams.uMax = width;
    stereoParams.vMax = height;
    stereoParams.setEqualMargin();
    
    
    ImageGenerator generator(&camera, foreImg, 250);
//    generator.setBackground(backImg);
    const int iterMax = root.get<int>("steps");
    int boardPoseCount = 0;
    const string imageBaseName = root.get<string>("output_name");
    
    results.open(root.get<string>("analysis_output_name"));
    
    for (auto & boardPoseItem : root.get_child("plane_transform"))
    {   
        int cameraIncCount = 0;
        generator.setPlaneTransform(readTransform(boardPoseItem.second));
        //depth GT
        Mat32f depthGT, depth, sigmaMat;
        generator.generateDepth(depthGT, xiCam0);
        
        imshow("depthGT", depthGT / 10);
        //base frame
        const string imgName = imageBaseName + "_" + to_string(boardPoseCount) + "_base.png";
        Mat8u img1 = imread(imgName, 0);
        Mat8u noise(img1.size());

        
        const int noiseLevel = root.get<int>("noise");
        if (noiseLevel != 0)
        {
            Mat8u noise(img1.size());
            randu(noise, 0, noiseLevel);
            img1 -= noise;
        }

        //different increment diretion
        for (auto & cameraIncItem : root.get_child("camera_increment"))
        {
            Transf dxi = readTransform(cameraIncItem.second);
            Transf xiCam = xiCam0.compose(dxi);
            cout << boardPoseCount << " " << cameraIncCount << endl;
            //increment count
            for (int i = 0; i < iterMax; i++, xiCam = xiCam.compose(dxi))
            {
                const string imgName = imageBaseName + "_" + to_string(boardPoseCount) 
                    + "_" + to_string(cameraIncCount) + "_" + to_string(i+1) + ".png";
                Mat8u img2 = imread(imgName, 0);

                if (noiseLevel != 0)
                {
                    Mat8u noise(img2.size());
                    randu(noise, 0, noiseLevel);
                    img2 -= noise;
                }
                Transf TleftRight = xiCam0.inverseCompose(xiCam);
                EnhancedSGM stereo(TleftRight, &camera, &camera, stereoParams);
                results << TleftRight.trans().norm() << "    ";
                DepthMap depthStereo;
                stereo.computeStereo(img1, img2, depthStereo);
                depthStereo.toMat(depth);
                depthStereo.sigmaToMat(sigmaMat);
                analyzeError(depthGT, depth, sigmaMat, stereoParams);
                results << endl;
                imshow("depth", depth / 10);
                imshow("img", img2);
                waitKey(30);
            }
            cameraIncCount++;
        }
        boardPoseCount++;
    }
    results.close();
    return 0;
}
