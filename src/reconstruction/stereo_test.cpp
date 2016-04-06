#include <iostream>
#include <fstream>
#include <ctime>

#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>

#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_stereo.h"

using namespace cv;
using namespace std;

typedef cv::Mat_<uint8_t> Mat8;

struct Ellipse
{
    double operator() (int x, int y)const
    {
        return (x - 50)*(x - 50) + 2*(y-50)*(y-50) - 2000;
    }
    
    double gradx(int x, int y)const
    {
        return 2*(x - 50);
    }
    
    double grady(int x, int y)const
    {
        return 4*(y - 50);
    }
};

int main()
{	
/*
    Matf test(100, 100);
    Matf test2;
    test.setTo(0);
    
    //(x - 20)^2 + 2(y-20)^2  = 100
    
    int x = 16;
    int y = 27;
    
    test(y, x) = 1;
    
    int x2 = 50;
    int y2 = 0;
    CurveRasterizer<Ellipse> curve(x, y, 0, 1, Ellipse());
    for (int i = 0; i < 3000; i++)
    {
        curve.makeStep();
        test *=0.99;
        test(curve.y, curve.x) = 1;
//        resize(test, test2, Size(0, 0), 4, 4);
        imshow("res", test);
        waitKey(35);
    }*/
     
    const double K = 1;
    array<double, 6> params1{0.571, 1.180, 378.304/K, 377.960/K, 654.923/K, 474.835/K};
    array<double, 6> params2{0.570, 1.186, 377.262/K, 376.938/K, 659.914/K, 489.024/K};
    Transformation<double> TleftRight(0.788019, 0.00459233, -0.0203431, -0.00243736, 0.0859855, 0.000375454);
    
//    Transformation<double> TleftRight( 0.5,    0.1, 0.0005, -0.025 , 0.005 ,  0.01);
//    array<double, 6> params1{0.5, 1, 150, 150, 320, 240};
//    array<double, 6> params2{0.5, 1, 150, 150, 320, 240};
//   
////    Transformation<double> TleftRight(1, 0, 0, 0, 0, 0);
//    Mat8 img1 = imread("/home/bogdan/projects/stack/perception/generate/res0.png", 0);
//    Mat8 img2 = imread("/home/bogdan/projects/stack/perception/generate/res5.png", 0);
    
    Mat8 img1 = imread("/home/bogdan/projects/data/icars/img_start/l01.jpg", 0);
    Mat8 img2 = imread("/home/bogdan/projects/data/icars/img_start/r01.jpg", 0);
    
//    GaussianBlur(img1, img1, Size(0, 0), 1, 1);
//    GaussianBlur(img2, img2, Size(0, 0), 1, 1);
//    EUCM cam1(params1.data());
//    EUCM cam2(params2.data());
    StereoEUCM stereo(TleftRight, img1.cols, img1.rows, params1.data(), params2.data());

    vector<Point> vec = {Point(932, 332), Point(1078, 348), Point(656, 534), Point(739,437)};
//    for (auto & pt : vec)
//    {
    cv::Mat_<uint8_t> out;
    
    auto t1 = clock();
    stereo.computeCost(img1, img2, out);
    auto t2 = clock();
    cv::Mat_<int> dynOut;
    
    cout << double(t2 - t1) / CLOCKS_PER_SEC << endl;
    
    stereo.computeDynamicProgramming(out, dynOut);
    auto t3 = clock();
    
    cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
    int dispMax = 64;
    int blockSize = 5;
    cv::Mat_<uint8_t> res(Size(out.cols/dispMax, out.rows));
    res.setTo(0);
    for (int v = 100; v < img1.rows - 100; v += blockSize)
    {
        auto dynRow = dynOut.row(v / blockSize);
        auto dynRow2 = dynOut.row(out.rows + v / blockSize);
        auto dynRow3 = dynOut.row(out.rows*2 + v / blockSize);
        auto dynRow4 = dynOut.row(out.rows*3 + v / blockSize);
//        auto dynRow2 = out.row(v / blockSize);
        auto resRow = res.row(v / blockSize);
        for (int u = 0; u < resRow.cols; u++)
        {
            unsigned int minVal = 1000000000;
            for (int i = 0; i < dispMax; i++)
            {
                unsigned int val = (dynRow(u*dispMax + i) + dynRow2(u*dispMax + i)) *0
                            + dynRow3(u*dispMax + i) + dynRow4(u*dispMax + i);
                if (minVal > val)
                {
                    resRow(u) = i*2;
                    minVal = val;
                    
                }
            }
//            cout << u << " " << v << " " << minVal << endl;
        }
         
    }
    resize(res, res, Size(0, 0), 5, 5, 0);
    imshow("res", res);
    waitKey(); 
//    for (int v = 100; v < img1.rows - 100; v++)
//    {
//        int idx = img1.cols * v + dispMax;
//        for (int u = dispMax; u < img1.cols - dispMax; u++, idx++)
//        {
//            int minVal = 200;
//            uint8_t * outPtr1 = out.data + (idx + 1)*dispMax;
//            uint8_t * outPtr2 = out.data + (idx - 1)*dispMax;
//            uint8_t * outPtr3 = out.data + (idx + img1.cols)*dispMax;
//            uint8_t * outPtr4 = out.data + (idx - img1.cols)*dispMax;
//            uint8_t * outPtr0 = out.data + idx*dispMax;
//            
//            uint8_t * resPtr = res.data + idx;
//            for (int i = 0; i < dispMax; i++, outPtr0++,
//                    outPtr1++, outPtr2++, outPtr3++, outPtr4++)
//            {
//                int err = *outPtr0 + *outPtr1 + *outPtr2 + *outPtr3 + *outPtr4;
//                if (err < minVal)
//                {
//                    *resPtr = i*6;
//                    minVal = err;
//                }
//            }
//        }
//    }
//    imshow("out", out);

//    }
    return 0;
}



