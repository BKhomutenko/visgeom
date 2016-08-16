#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"
#include "reconstruction/eucm_motion_stereo.h"

int main(int argc, char** argv)
{	
    /*Polynomial2 poly2;
    poly2.kuu = -1; 
    poly2.kuv = 1; 
    poly2.kvv= -1; 
    poly2.ku = 0.25; 
    poly2.kv = 0.25; 
    poly2.k1 = 5;
    
    CurveRasterizer<Polynomial2> raster(1, 1, -100, 100, poly2);
    CurveRasterizer2<Polynomial2> raster2(1, 1, -100, 100, poly2);

    auto tr0 = clock();
    int x1 = 0;
    int x2 = 0;
    for (int i = 0; i < 10000000; i++)
    {
        raster.step();
        x1 += raster.x;
    }
    auto tr1 = clock();
    
    for (int i = 0; i < 10000000; i++)
    {
        raster2.step();
        x2 += raster2.x;
    }
    auto tr2 = clock();
    
    cout << "optimized " << double(tr1 - tr0) / CLOCKS_PER_SEC << endl;
    cout << "simple " << double(tr2 - tr1) / CLOCKS_PER_SEC << endl;
    cout << x1 << " " << x2 << endl;
    return 0;*/
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
    
    array<double, 6> robotPose3, robotPose4;
    cout << "Third robot's pose :" << endl;
    for (auto & e: robotPose3) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    cout << "Fourth robot's pose :" << endl;
    for (auto & e: robotPose4) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    Transformation<double> T03(robotPose3.data()), T04(robotPose4.data());
    Transformation<double> TleftRight1 = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    Transformation<double> TleftRight3 = T03.compose(TbaseCamera).inverseCompose(T04.compose(TbaseCamera));
    StereoParameters stereoParams;
    stereoParams.verbosity = 1;
    paramFile >> stereoParams.u0;
    paramFile >> stereoParams.v0;
    paramFile >> stereoParams.dispMax;
    paramFile >> stereoParams.scale;
    paramFile.ignore();
//    stereoParams.lambdaStep = 3;
//    stereoParams.lambdaJump = 15;
    string fileName1, fileName2;
    string fileName3, fileName4;
    getline(paramFile, fileName1);
    getline(paramFile, fileName2);
    getline(paramFile, fileName3);
    getline(paramFile, fileName4);
    
    Mat8u img1 = imread(fileName1, 0);
    Mat8u img2 = imread(fileName2, 0);
    Mat8u img3 = imread(fileName1, 0);
    Mat8u img4 = imread(fileName2, 0);
    
    stereoParams.uMax = img1.cols;
    stereoParams.vMax = img1.rows;
    stereoParams.setEqualMargin();
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
    EnhancedStereo stereo1(TleftRight1, &camera1, &camera2, stereoParams);
    EnhancedStereo stereo3(TleftRight3, &camera1, &camera2, stereoParams);
    cout << "    initialization time : " << timer.elapsed() << endl;
    
    DepthMap depth1, depth3;
    
    MotionStereoParameters motionParams;
    
    MotionStereo motionStereo(&camera1, &camera2, motionParams);
    
    timer.reset();
    
    stereo1.computeStereo(img1, img2, depth1);
    
    cout << "    first stereo : " << timer.elapsed() << endl;
    
    Mat32f depth1Mat, depth2Mat;
    depth1.toMat(depth1Mat);
    
    timer.reset();
    
    motionStereo.setBaseImage(img1);
    Transformation<double> T13 = T01.compose(TbaseCamera).inverseCompose(T03.compose(TbaseCamera));
    motionStereo.reprojectDepth(T13, img3, depth1);
    cout << "    reproject depth : " << timer.elapsed() << endl;
    
    depth1.toMat(depth2Mat);
    
    
    //TODO compare to the ground truth
//    
//    for (auto & x : {Point(320, 300), Point(500, 300), Point(750, 300), Point(350, 500), Point(600, 450)})
//    {
//        out1(x) = 0;
//        out1(x.y + 1, x.x) = 0;
//        out1(x.y, x.x + 1) = 255;
//        out1(x.y + 1, x.x + 1) = 255;
//        stereo.traceEpipolarLine(x, out2);
//    }
    
    
   
    imshow("out1", depth1Mat / 3);
    imshow("out2", depth2Mat / 3);
    waitKey(); 
    return 0;
}



