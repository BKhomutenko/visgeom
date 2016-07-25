#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "reconstruction/curve_rasterizer.h"
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
    paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    
    Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    
    MotionStereoParameters stereoParams;
    stereoParams.verbosity = 1;
    int foo;
    paramFile >> foo;
    paramFile >> foo;
    paramFile >> foo;
    paramFile >> stereoParams.scale;
    paramFile.ignore();
    string fileName1, fileName2;
    getline(paramFile, fileName1);
    getline(paramFile, fileName2);
    
    Mat8u img1 = imread(fileName1, 0);
    Mat8u img2 = imread(fileName2, 0);
    Mat16s img1lap, img2lap;

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
    stereoParams.verbosity = 0;
    MotionStereo stereo(&camera1, &camera2, stereoParams);
    stereo.setBaseImage(img1);
    cout << "    initialization time : " << timer.elapsed() << endl;
    Mat8u out1, out2;
    
    img1.copyTo(out1);
    img2.copyTo(out2);
    
//    
//    for (auto & x : {Point(320, 300), Point(500, 300), Point(750, 300), Point(350, 500), Point(600, 450)})
//    {
//        out1(x) = 0;
//        out1(x.y + 1, x.x) = 0;
//        out1(x.y, x.x + 1) = 255;
//        out1(x.y + 1, x.x + 1) = 255;
//        stereo.traceEpipolarLine(x, out2);
//    }
    
    DepthMap depthOut;
    
    int u0 = 25;
    int v0 = 25;
    int w = (img1.cols - 2*u0) / stereoParams.scale;
    int h = (img1.rows - 2*u0) / stereoParams.scale;
    
    
    DepthMap depth(&camera1, w, h, u0, v0, stereoParams.scale);
    depth.setDefault();
    Mat8u res;
    timer.reset();
    stereo.computeDepth(TleftRight, img2, depth);
    cout << timer.elapsed() << endl;
    timer.reset();

    res.create(h, w);
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            res(y, x) = depth.at(x, y);
        }
    }
    
    
    imshow("out1", out1);
    imshow("out2", out2);
    imshow("res", res*100);
    waitKey(); 
    return 0;
}



