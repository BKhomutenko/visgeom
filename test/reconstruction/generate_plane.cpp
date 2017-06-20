#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "projection/eucm.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

using namespace std;
using namespace cv;

//void showHistogram(const string & name, const Mat & img)
//{
//    int histSize = 256;    // bin size
//    float range[] = { 0, 1 };
//    const float *ranges[] = { range };

//    // Calculate histogram
//    MatND hist;
//    calcHist( &img, 1, 0, Mat(), hist, 1, &histSize, ranges, true, false );
//    
//    // Show the calculated histogram in command window
//    double total;
//    total = img.rows * img.cols;
//    for( int h = 0; h < histSize; h++ )
//         {
//            float binVal = hist.at<float>(h);
//            cout<<" "<<binVal;
//         }

//    // Plot the histogram
//    int hist_w = 512; int hist_h = 400;
//    int bin_w = round( (double) hist_w/histSize );

//    Mat histImage( hist_h, hist_w, CV_8UC1, Scalar( 0,0,0) );
//    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
//    
//    for( int i = 1; i < histSize; i++ )
//    {
//      line( histImage, Point( bin_w*(i-1), hist_h - round(hist.at<float>(i-1)) ) ,
//                       Point( bin_w*(i), hist_h - round(hist.at<float>(i)) ),
//                       Scalar( 255, 0, 0), 2, 8, 0  );
//    }

//    imshow(name, histImage );
//}

int main(int argc, char** argv)
{	
    
    ScaleParameters params;
    params.scale = 2;
    params.u0 = 0; params.v0 = 0;
    params.uMax = 600; params.vMax = 400; //image size
    params.xMax = 300; params.yMax = 200; //scaled size
    vector<double> cameraParams{0.5, 1, 250, 250, 300, 200};
    EnhancedCamera camera(cameraParams.data());
    const double L = 0.3;
    Vector3dVec polygonVec{ Vector3d(   -L,     -L,     0),
                            Vector3d(    L,     -L,     0),
                            Vector3d(    L,      L,     0),
                            Vector3d(   -L,      L,     0) };
    DepthMap depth = DepthMap::generatePlane(&camera, params, Transf(0, 0, 0.5, 0.1, 0.1, 0.1), polygonVec);
    
    Mat32f depthMat, depthMat2;
    depth.toMat(depthMat);
    DepthMap depth2 = depth.wrapDepth(Transf(0, 0, 0, 0.1, 0.1, 0.1));
    depth2.toMat(depthMat2);
    imshow("depth", depthMat);
    imshow("depth2", depthMat2);           
    waitKey();
    return 0;
}



