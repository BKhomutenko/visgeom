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

#include "render/render.h"
#include "render/background.h"
#include "render/plane.h"

//FIXME make an argument

int main(int argc, char** argv) 
{
//    ptree root;
//    read_json(argv[1], root);
//    const vector<double> intrinsic = readVector(root.get_child("camera_intrinsics"));
//    const int width = root.get<int>("image.width");
//    const int height = root.get<int>("image.height");
//    Transf xiCam0 = readTransform(root.get_child("camera_transform"));
//    
//    Mat8u foreImg = imread(root.get<string>("foreground"), 0);
//    
//    EnhancedCamera camera(width, height, intrinsic.data());
//    
//    //init stereoParameters
//    SGMParameters stereoParams;
//    
//    stereoParams.verbosity = root.get<int>("stereo.verbosity");
//    stereoParams.salientPoints = true;
//    stereoParams.u0 = root.get<int>("stereo.u0");
//    stereoParams.v0 = root.get<int>("stereo.v0");
//    stereoParams.dispMax = root.get<int>("stereo.disparity_max");
//    stereoParams.scale = root.get<int>("stereo.scale");
//    stereoParams.flawCost = root.get<int>("stereo.flaw_cost");
//    stereoParams.scaleVec = readIntVector(root.get_child("stereo.scale_vector"));
//    stereoParams.uMax = width;
//    stereoParams.vMax = height;
//    stereoParams.setEqualMargin();
//    
//    
//    ImageGenerator generator(&camera, foreImg, 250);
////    generator.setBackground(backImg);
//    const int iterMax = root.get<int>("steps");
//    int boardPoseCount = 0;
//    const string imageBaseName = root.get<string>("output_name");
//    
//    results.open(root.get<string>("analysis_output_name"));
//    
//    for (auto & boardPoseItem : root.get_child("plane_transform"))
//    {   
//        int cameraIncCount = 0;
//        generator.setPlaneTransform(readTransform(boardPoseItem.second));
//        //depth GT
//        Mat32f depthGT, depth, sigmaMat;
//        generator.generateDepth(depthGT, xiCam0);
//        
//        imshow("depthGT", depthGT / 10);
//        //base frame
//        const string imgName = imageBaseName + "_" + to_string(boardPoseCount) + "_base.png";
//        Mat8u img1 = imread(imgName, 0);
//        Mat8u noise(img1.size());
//        
//        const int noiseLevel = root.get<int>("noise");
//        if (noiseLevel != 0)
//        {
//            Mat8u noise(img1.size());
//            randu(noise, 0, noiseLevel);
//            img1 -= noise;
//        }
//        //different increment diretion
//        for (auto & cameraIncItem : root.get_child("camera_increment"))
//        {
//            Transf dxi = readTransform(cameraIncItem.second);
//            Transf xiCam = xiCam0.compose(dxi);
//            cout << boardPoseCount << " " << cameraIncCount << endl;
//            //increment count
//            for (int i = 0; i < iterMax; i++, xiCam = xiCam.compose(dxi))
//            {
//                const string imgName = imageBaseName + "_" + to_string(boardPoseCount) 
//                    + "_" + to_string(cameraIncCount) + "_" + to_string(i+1) + ".png";
//                Mat8u img2 = imread(imgName, 0);
//                if (noiseLevel != 0)
//                {
//                    Mat8u noise(img2.size());
//                    randu(noise, 0, noiseLevel);
//                    img2 -= noise;
//                }
//                Transf TleftRight = xiCam0.inverseCompose(xiCam);
//                EnhancedSGM stereo(TleftRight, &camera, &camera, stereoParams);
//                results << TleftRight.trans().norm() << "    ";
//                DepthMap depthStereo;
//                stereo.computeStereo(img1, img2, depthStereo);
//                depthStereo.toMat(depth);
//                depthStereo.sigmaToMat(sigmaMat);
//                analyzeError(depthGT, depth, sigmaMat, stereoParams);
//                results << endl;
//                imshow("depth", depth / 10);
//                imshow("img", img2);
//                waitKey(30);
//            }
//            cameraIncCount++;
//        }
//        boardPoseCount++;
//    }
//    results.close();
    array<double, 6> params {0.5, 1, 150, 150, 320, 240};
    EnhancedCamera camera(640, 480, params.data());
    Renderer renderer(&camera);
    
    Mat8u img = imread("/home/bogdan/Pictures/space-wallpaper.png", 0);
    
    Background * back = new Background(img);
    back->_u0 = img.cols / 2;
    back->_v0 = img.rows / 2;
    back->_fu = 200;
    back->_fv = 200;
    
    renderer._objectVec.push_back(back);
    
    img = imread("/home/bogdan/projects/data/render2/foreground.jpg", 0);
    Plane * plane = new Plane(Transf(0, 0.6, 0, -M_PI/2, 0, 0), img);
    plane->_u0 = img.cols / 2;
    plane->_v0 = img.rows / 2;
    plane->_fu = 100;
    plane->_fv = 100;
    renderer._objectVec.push_back(plane);
    
    img = imread("/home/bogdan/projects/stack/aliasing.png", 0);
//    img = imread("/home/bogdan/projects/data/render/background.jpg", 0);
    plane = new Plane(Transf(0, 0, 1, 0, 0, 0), img);
    plane->_u0 = img.cols / 2;
    plane->_v0 = img.rows / 2;
    plane->_fu = 100;
    plane->_fv = 100;
    renderer._objectVec.push_back(plane);
    
    renderer.fillBuffers();
    Mat8u res;
    renderer.fillImage(res);
    imshow("depth", renderer._depthMat/5);
    imshow("res", res);
    waitKey();
    return 0;
}
