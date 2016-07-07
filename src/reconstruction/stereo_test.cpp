#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"

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
    paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    
    Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    
    StereoParameters stereoParams;
    paramFile >> stereoParams.uMargin;
    paramFile >> stereoParams.vMargin;
    paramFile >> stereoParams.dispMax;
    paramFile >> stereoParams.scale;
    paramFile.ignore();
//    stereoParams.lambdaStep = 5;
//    stereoParams.lambdaJump = 15;
    string fileName1, fileName2;
    getline(paramFile, fileName1);
    getline(paramFile, fileName2);
    
    Mat8u img1 = imread(fileName1, 0);
    Mat8u img2 = imread(fileName2, 0);
    Mat16s img1lap, img2lap;

    stereoParams.imageWidth = img1.cols;
    stereoParams.imageHeight = img1.rows;
//    
//    Laplacian(img1, img1lap, CV_16S, 3);
//    Laplacian(img2, img2lap, CV_16S, 3);
    
//    GaussianBlur(img1, img1, Size(5, 5), 0, 0);
//    GaussianBlur(img2, img2, Size(5, 5), 0, 0);
    
//    img1lap = img1lap + 128;
//    img2lap = img2lap + 128;
//    img1lap.copyTo(img1);
//    img2lap.copyTo(img2);
//    
    EnhancedStereo stereo(TleftRight, params1.data(), params2.data(), stereoParams);
    
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
    

    Mat8u res;
    auto t2 = clock();
    stereo.computeCost(img1, img2);
    auto t3 = clock();
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    t2 = clock();
    stereo.computeDynamicProgramming();
    t3 = clock();
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    t2 = clock();
    stereo.reconstructDisparity();
    t3 = clock();
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    t2 = clock();
    stereo.upsampleDisparity(img1, res);
    t3 = clock();
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    
    imshow("out1", out1);
    imshow("out2", out2);
    imshow("res", res);
    waitKey(); 
    return 0;
}



