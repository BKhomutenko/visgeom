#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/depth_map.h"

EnhancedEpipolar * epipolar;
Mat8u img1, img2;

void CallBackFunc1(int event, int x, int y, int flags, void* userdata)
{
     if  ( event == cv::EVENT_LBUTTONDOWN )
     {
          cout << "Left button of the mouse is clicked 1 - position (" << x << ", " << y << ")" << endl;
          epipolar->traceEpipolarLine(x, y, img2, CAMERA_1, 100);
          cv::circle(img1, Point(x, y), 0, Scalar(0), -1);
          imshow("out1", img1);
          imshow("out2", img2);
     }
}

void CallBackFunc2(int event, int x, int y, int flags, void* userdata)
{
     if  ( event == cv::EVENT_LBUTTONDOWN )
     {
          cout << "Left button of the mouse is clicked 2 - position (" << x << ", " << y << ")" << endl;
          epipolar->traceEpipolarLine(x, y, img1,  CAMERA_2, 100);
          cv::circle(img2, Point(x, y), 0, Scalar(0), -1);
          imshow("out1", img1);
          imshow("out2", img2);
     }
}

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
    
    EnhancedCamera cam1(params1.data()), cam2(params2.data());
    
    epipolar = new EnhancedEpipolar(TleftRight, &cam1, &cam2, 2000);
    
    string fileName1, fileName2;
    
    getline(paramFile, fileName1); //to SGM parameters
    
    while(getline(paramFile, fileName1))
    {
        getline(paramFile, fileName2);
        
        img1 = imread(fileName1, 0);
        img2 = imread(fileName2, 0);
        
        imshow("out1", img1);
        imshow("out2", img2);
        cv::setMouseCallback("out1", CallBackFunc1, NULL);
        cv::setMouseCallback("out2", CallBackFunc2, NULL);
        waitKey(); 
    }
    delete epipolar;
    return 0;
}



