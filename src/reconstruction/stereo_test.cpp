#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>

#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"

using namespace cv;
using namespace std;

typedef cv::Mat_<uint8_t> Mat8;

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
    
    array<double, 6> transfArr;
    cout << "Transformation First-Second:" << endl;
    for (auto & e: transfArr) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TleftRight(transfArr.data());
    
    StereoParameters stereoParams;
    paramFile >> stereoParams.u0;
    paramFile >> stereoParams.v0;
    paramFile >> stereoParams.disparityMax;
    paramFile >> stereoParams.blockSize;
    paramFile.ignore();
    
    string fileName1, fileName2;
    getline(paramFile, fileName1);
    getline(paramFile, fileName2);
    
    Mat8 img1 = imread(fileName1, 0);
    Mat8 img2 = imread(fileName2, 0);
  
//    const double K = 1;
//    array<double, 6> params1{0.571, 1.180, 378.304/K, 377.960/K, 654.923/K, 474.835/K};
//    array<double, 6> params2{0.570, 1.186, 377.262/K, 376.938/K, 659.914/K, 489.024/K};
//    Transformation<double> TleftRight(0.788019, 0.00459233, -0.0203431, -0.00243736, 0.0859855, 0.000375454);
//    Mat8 img1 = imread("/home/bogdan/projects/data/icars/img_left/frame0010.jpg", 0);
//    Mat8 img2 = imread("/home/bogdan/projects/data/icars/img_right/frame0010.jpg", 0);
//    Mat8 img1 = imread("/home/bogdan/projects/data/icars/img_start/l01.jpg", 0);
//    Mat8 img2 = imread("/home/bogdan/projects/data/icars/img_start/r01.jpg", 0);    
    
//    array<double, 6> params1{0.01, 1, 500, 500, 1800/4, 1500/4};
//    array<double, 6> params2{0.01, 1, 500, 500, 1800/4, 1500/4};
//    Transformation<double> TleftRight(1, 0, 0, 0, 0, 0);
//    
//    resize(img1, img1, Size(0, 0), 0.5, 0.5);
//    resize(img2, img2, Size(0, 0), 0.5, 0.5);

    
//    Transformation<double> TleftRight( 0.5,    0.1, 0.0005, -0.025 , 0.005 ,  0.01);
//    array<double, 6> params1{0.5, 1, 150, 150, 320, 240};
//    array<double, 6> params2{0.5, 1, 150, 150, 320, 240};
//    Mat8 img1 = imread("/home/bogdan/projects/stack/perception/generate/res0.png", 0);
//    Mat8 img2 = imread("/home/bogdan/projects/stack/perception/generate/res5.png", 0);
    
    
    
//    GaussianBlur(img1, img1, Size(0, 0), 1, 1);
//    GaussianBlur(img2, img2, Size(0, 0), 1, 1);
//    EUCM cam1(params1.data());
//    EUCM cam2(params2.data());
    
    EnhancedStereo stereo(TleftRight, img1.cols, img1.rows, params1.data(), params2.data(), stereoParams);

    cv::Mat_<uint8_t> res;
    auto t2 = clock();
    stereo.comuteStereo(img1, img2, res);
    auto t3 = clock();
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    
    imshow("res", res);
    waitKey(); 
    return 0;
}



