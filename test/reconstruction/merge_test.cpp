#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"

int main(int argc, char** argv)
{	

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
    SGMParameters stereoParams;
    stereoParams.verbosity = 1;
    stereoParams.hypMax = 1;
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
    Mat8u img3 = imread(fileName3, 0);
    Mat8u img4 = imread(fileName4, 0);
    
    stereoParams.uMax = img1.cols;
    stereoParams.vMax = img1.rows;
    stereoParams.setEqualMargin();
    
    Timer timer, timer1;
    timer1.reset();
    EnhancedCamera camera1(params1.data()), camera2(params2.data());
    EnhancedSGM stereo1(TleftRight1, &camera1, &camera2, stereoParams);
    EnhancedSGM stereo3(TleftRight3, &camera1, &camera2, stereoParams);
    cout << "    initialization time : " << timer.elapsed() << endl;
    
    DepthMap depth1, depth3;
    
    MotionStereoParameters motionParams;
    motionParams.dispMax = 96;
    MotionStereo motionStereo(&camera1, &camera2, motionParams);
    
    timer.reset();
    
    stereo1.computeStereo(img1, img2, depth1);
    
    cout << "    first stereo : " << timer.elapsed() << endl;
    
    Mat32f depth1Mat, depth2Mat, depth3Mat, depth3FilteredMat;
    depth1.toMat(depth1Mat);
    
    timer.reset();
    
    motionStereo.setBaseImage(img1);
    Transformation<double> T13 = T01.compose(TbaseCamera).inverseCompose(T03.compose(TbaseCamera));
    motionStereo.reprojectDepth(T13, img3, depth1);
    cout << "    reproject depth : " << timer.elapsed() << endl;
    
    depth1.toMat(depth2Mat);
    
    timer.reset();
    
    stereo1.computeStereo(img3, img4, depth3);
    
    cout << "    second stereo : " << timer.elapsed() << endl;
    
    timer.reset();
    depth3.merge(depth1);
    cout << "    merge : " << timer.elapsed() << endl;
    depth3.toMat(depth3Mat);

//    depth3.filterNoise();
//    depth3.toMat(depth3FilteredMat);


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
    
    
   
    
    imshow("img1", img1);
    imshow("img3", img3);
    imshow("out1", depth1Mat / 3);
    imshow("out2", depth2Mat / 3);
    imshow("out3", depth3Mat / 3);
//    imshow("out3", depth3FilteredMat / 3);

    waitKey(); 
    return 0;
}



