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
#include "localization/photometric.h"

#include "io.h"
#include "ocv.h"
#include "timer.h"

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/depth_map.h"

using namespace cv;
using namespace std;

// TODO separate the ground truth generation from the stereo 
int main (int argc, char const* argv[])
{

/// RECONSTRUCTION PART
    ifstream paramFile(argv[1]);
    if (not paramFile.is_open())
    {
        cout << argv[1] << " : ERROR, file is not found" << endl;
        return 0;
    }
    
    array<double, 6> params1;
    array<double, 6> params2;
    
    cout << "First EU Camera model parameters :" << endl;
    for (auto & p: params1) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    cout << "Second EU Camera model parameters :" << endl;
    for (auto & p: params2) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 6> cameraPose;
    cout << "Camera pose wrt the robot :" << endl;
    for (auto & e: cameraPose) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbaseCamera(cameraPose.data());
    
    array<double, 6> robotPose1, robotPose2;
    cout << "First robot's pose :" << endl;
    for (auto & e: robotPose1) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    cout << "Second robot's pose :" << endl;
    for (auto & e: robotPose2) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    
    Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    
    SGMParameters stereoParams;
    stereoParams.flawCost = 5;
    stereoParams.verbosity = 1;
    stereoParams.hypMax = 3;
    stereoParams.salientPoints = true;
    paramFile >> stereoParams.u0;
    paramFile >> stereoParams.v0;
    paramFile >> stereoParams.dispMax;
    paramFile >> stereoParams.scale;
    paramFile.ignore();
//    stereoParams.lambdaStep = 3;
//    stereoParams.lambdaJump = 15;
    string fileName1, fileName2;
    getline(paramFile, fileName1);
    getline(paramFile, fileName2);
    
    Mat8u img1 = imread(fileName1, 0);
    Mat8u img2 = imread(fileName2, 0);
    Mat16s img1lap, img2lap;

    stereoParams.uMax = img1.cols;
    stereoParams.vMax = img1.rows;
    stereoParams.setXMargin(stereoParams.u0);
    stereoParams.setYMargin(300);
//    
//    Laplacian(img1, img1lap, CV_16S, 3);
//    Laplacian(img2, img2lap, CV_16S, 3);
//    
//    GaussianBlur(img1, img1, Size(3, 3), 0, 0);
//    GaussianBlur(img2, img2, Size(3, 3), 0, 0);
//    
//    img1lap = img1lap + 128;
//    img2lap = img2lap + 128;
//    img1lap.copyTo(img1);
//    img2lap.copyTo(img2);
////    
    
    Timer timer;
    EnhancedCamera camera1(params1.data()), camera2(params2.data());
    EnhancedSGM stereo(TleftRight, &camera1, &camera2, stereoParams);
    cout << "    initialization time : " << timer.elapsed() << endl;
    
//    
//    for (auto & x : {Point(320, 300), Point(500, 300), Point(750, 300), Point(350, 500), Point(600, 450)})
//    {
//        out1(x) = 0;
//        out1(x.y + 1, x.x) = 0;
//        out1(x.y, x.x + 1) = 255;
//        out1(x.y + 1, x.x + 1) = 255;
//        stereo.traceEpipolarLine(x, out2);
//    }
    

//    Mat32f res;
//    timer.reset();
//    stereo.computeCurveCost(img1, img2);
//    cout << timer.elapsed() << endl;
//    timer.reset();
//    stereo.computeDynamicProgramming();
//    cout << timer.elapsed() << endl;
//    timer.reset();
//    stereo.reconstructDisparity();
//    cout << timer.elapsed() << endl;
//    timer.reset();
//    stereo.computeDepth(res);
//    cout << timer.elapsed() << endl;
    
    DepthMap depth;
    timer.reset();
    
//    stereo.computeStereo(img1, img2, depthMat);
    stereo.computeStereo(img1, img2, depth);
    cout << "    stereo total time : " << timer.elapsed() << endl;
    
    
    
    paramFile.close();
    
/// LOCALIZATION PART
    paramFile.open(argv[2]);
    if (not paramFile.is_open())
    {
        cout << argv[2] << " : ERROR, file is not found" << endl;
        return 0;
    }
    
    array<double, 6> params;
    
    cout << "EU Camera model parameters :" << endl;
    for (auto & p: params) 
    {
        paramFile >> p;
        cout << setw(15) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 6> cameraPose0;
    cout << "Camera pose wrt the robot :" << endl;
    for (auto & e: cameraPose0) 
    {
        paramFile >> e;
        cout << setw(15) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbaseCamera1(cameraPose0.data());
    
    double foo;
    ScaleParameters scaleParams;
    paramFile >> scaleParams.u0;
    paramFile >> scaleParams.v0;
    paramFile >> foo;
    paramFile >> scaleParams.scale;
    paramFile.ignore();
    
    string imageDir;
    getline(paramFile, imageDir);
    
    //read the first image name and the robot's pose
    string imageInfo, imageName;
    getline(paramFile, imageInfo);
    istringstream imageStream(imageInfo);
    imageStream >> imageName;
    for (auto & x : robotPose1) imageStream >> x;
    
    Mat32f img11 = imread(imageDir + imageName, 0);
    
    scaleParams.uMax = img11.cols;
    scaleParams.vMax = img11.rows;
    scaleParams.setEqualMargin();
    
    // create a camera
    EnhancedCamera camera(params.data());
    
    // Init the distance map
    T01 = Transformation<double>(robotPose1.data());
    Transformation<double> T0Camera = T01.compose(TbaseCamera1);
    
//     Init the localizer
    ScalePhotometric localizer(5, &camera);
    localizer.setVerbosity(1);
    localizer.computeBaseScaleSpace(img11);
    localizer.depth() = depth;
    
     
    imshow("img1", img11/255);
    
    std::chrono::high_resolution_clock clock_;
    
    while (getline(paramFile, imageInfo))
    {
        istringstream imageStream(imageInfo);
        
        imageStream >> imageName;
        for (auto & x : robotPose2) imageStream >> x;
    
        Transformation<double> T02(robotPose2.data());
        Transformation<double> T12 = T0Camera.inverseCompose(T02.compose(TbaseCamera1));
        cout << " WO POSE : " << T12 << endl;
//        T12 = T12.compose(Transformation<double>(-0.01, -0.01, -0.3, -0.003, -0.003, -0.005));
//        T12 = T12.compose(Transformation<double>(-0.001, 0.005, -0.01, 0.001, 0.001, 0.001));
        Mat32f img2 = imread(imageDir + imageName, 0);
        cout << " img2 : " << imageName << endl;
        imshow("img2", img2/255);
        Mat32f img2WrapOrig, img2WrapFinal;
        
        localizer.wrapImage(img2, img2WrapOrig, T12);
        for (int iter = 0; iter < 1; iter++)
        {
            Timer timer;
            localizer.computePose(img2, T12);
            cout << timer.elapsed() << endl;
            cout << T12 << endl;
            
            localizer.wrapImage(img2, img2WrapFinal, T12);
            imshow("img2WrapOrig", img2WrapOrig/256);
            imshow("img2WrapFinal", img2WrapFinal/256);
            
            imshow("delta img2WrapOrig", abs(img11 - img2WrapOrig)/250);
            imshow("delta img2WrapFinal", abs(img11 - img2WrapFinal)/250);
            waitKey();
        }
    }
    return 0;
}
