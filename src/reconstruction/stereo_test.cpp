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



