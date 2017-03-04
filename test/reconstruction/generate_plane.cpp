#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "projection/eucm.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

using namespace std;
using namespace cv;

void showHistogram(const string & name, const Mat & img)
{
    int histSize = 256;    // bin size
    float range[] = { 0, 1 };
    const float *ranges[] = { range };

    // Calculate histogram
    MatND hist;
    calcHist( &img, 1, 0, Mat(), hist, 1, &histSize, ranges, true, false );
    
    // Show the calculated histogram in command window
    double total;
    total = img.rows * img.cols;
    for( int h = 0; h < histSize; h++ )
         {
            float binVal = hist.at<float>(h);
            cout<<" "<<binVal;
         }

    // Plot the histogram
    int hist_w = 512; int hist_h = 400;
    int bin_w = round( (double) hist_w/histSize );

    Mat histImage( hist_h, hist_w, CV_8UC1, Scalar( 0,0,0) );
    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    
    for( int i = 1; i < histSize; i++ )
    {
      line( histImage, Point( bin_w*(i-1), hist_h - round(hist.at<float>(i-1)) ) ,
                       Point( bin_w*(i), hist_h - round(hist.at<float>(i)) ),
                       Scalar( 255, 0, 0), 2, 8, 0  );
    }

    imshow(name, histImage );
}

int main(int argc, char** argv)
{	

    ifstream paramFile(argv[1]);
    if (not paramFile.is_open())
    {
        cout << argv[1] << " : ERROR, file is not found" << endl;
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
    
    array<double, 6> cameraPose;
    cout << "Camera pose wrt the robot :" << endl;
    for (auto & e: cameraPose) 
    {
        paramFile >> e;
        cout << setw(15) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbaseCamera(cameraPose.data());
    
    array<double, 6> planePose;
    cout << "Plane pose :" << endl;
    for (auto & e: cameraPose) 
    {
        paramFile >> e;
        cout << setw(15) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbasePlane(cameraPose.data());
    
    SGMParameters stereoParams;
    stereoParams.verbosity = 2;
    stereoParams.salientPoints = false;
    stereoParams.maxBias = 8;
    paramFile >> stereoParams.u0;
    paramFile >> stereoParams.v0;
    paramFile >> stereoParams.dispMax;
    paramFile >> stereoParams.scale;
    
    paramFile.ignore();
    
    string imageDir;
    getline(paramFile, imageDir);
    
    string imageInfo, imageName;
    array<double, 6> robotPose1, robotPose2;
    getline(paramFile, imageInfo);
    istringstream imageStream(imageInfo);
    imageStream >> imageName;
    for (auto & x : robotPose1) imageStream >> x;

    Mat8u img1 = imread(imageDir + imageName, 0);
    stereoParams.uMax = img1.cols;
    stereoParams.vMax = img1.rows;
    stereoParams.setEqualMargin();
    
    EnhancedCamera camera(params.data());
    int counter = 2;
    while (getline(paramFile, imageInfo))
    {
        istringstream imageStream(imageInfo);
        
        imageStream >> imageName;
        for (auto & x : robotPose2) imageStream >> x;
    
        Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
        Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
        
        Mat8u img2 = imread(imageDir + imageName, 0);
        EnhancedSGM stereo(TleftRight, &camera, &camera, stereoParams);

        DepthMap depthStereo;
        Mat32f distMat;
        auto t2 = clock();
        stereo.computeStereo(img1, img2, depthStereo);
        auto t3 = clock();
        
        depthStereo.toMat(distMat);
        
//        cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
        
        Transformation<double> T0Camera = T01.compose(TbaseCamera);
        DepthMap depth = DepthMap::generatePlane(&camera, stereoParams,
                T0Camera.inverseCompose(TbasePlane),
                Vector3dVec{Vector3d(-0.1, -0.1, 0), Vector3d(-0.1 + 3 * 0.45, -0.1, 0),
                          Vector3d(-0.1 + 3 * 0.45, 0.5, 0), Vector3d(-0.1, 0.5, 0) } );
        Mat32f planeMat;
        depth.toMat(planeMat);
        imshow("dist" + to_string(counter) , distMat/2);
        imshow("plane" , planeMat/2);
        imwrite("/home/bogdan/projects/plane.png", planeMat);
        double err = 0;
        double err2 = 0;
        double dist = 0;
        int N = 0;
        int Nmax = 0;
        Mat32f inlierMat(planeMat.size());
        inlierMat.setTo(0);
        for (int u = 0; u < distMat.cols; u++)
        {
            for (int v = 0; v < distMat.rows; v++)
            {
                
                if (planeMat(v, u) == 0) continue;
                Nmax++;
                dist += planeMat(v, u);
                inlierMat(v, u) = 1;
                if (distMat(v, u) == 0 or distMat(v, u) != distMat(v, u) or planeMat(v, u) != planeMat(v, u)) continue;
                if (abs(distMat(v, u) - planeMat(v, u)) > 0.10) continue;
                inlierMat(v, u) = 0;
                err += distMat(v, u) - planeMat(v, u);
                err2 += pow(distMat(v, u) - planeMat(v, u), 2);
                N++;
            }
        }
//        cout << (counter - 1) * 7 << " & " << dist/ Nmax * 1000 << " & " << err / N *1000 << " & " << sqrt(err2 / N)*1000  
//                << " & " << 100 * N / double(Nmax) << "\\\\" << endl << "\\hline" << endl;
        cout << "avg err : " << err / N *1000 << " avg err2 : " 
<< sqrt(err2 / N)*1000  << " number of inliers : " << 100 * N / double(Nmax) << endl;
        cout << planeMat.size() << " " <<  distMat.size() << endl;
        imshow("diff" + to_string(counter), abs(planeMat - distMat));
        imshow("inliers" + to_string(counter), inlierMat);
        counter++;
        
    }
    waitKey();
    return 0;
}



