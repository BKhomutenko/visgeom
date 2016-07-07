/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#include <iostream>
#include <iomanip>
#include <fstream>

#include <opencv2/opencv.hpp>

#include "localization/photometric.h"

using namespace cv;
using namespace std;

typedef cv::Mat_<float> Matf;

// TODO separate the ground truth generation from the stereo 
int main (int argc, char const* argv[])
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
    
    StereoParameters stereoParams;
    paramFile >> stereoParams.uMargin;
    paramFile >> stereoParams.vMargin;
    paramFile >> stereoParams.dispMax;
    paramFile >> stereoParams.scale;
    paramFile.ignore();
    
    string imageDir;
    getline(paramFile, imageDir);
    
    //read the first image name and the robot's pose
    string imageInfo, imageName;
    array<double, 6> robotPose1, robotPose2;
    getline(paramFile, imageInfo);
    istringstream imageStream(imageInfo);
    imageStream >> imageName;
    for (auto & x : robotPose1) imageStream >> x;
    
    Matf img1 = imread(imageDir + imageName, 0);
    
    stereoParams.imageWidth = img1.cols;
    stereoParams.imageHeight = img1.rows;
    
    const int sigma = 20;
    GaussianBlur(img1, img1, Size(0, 0), sigma, sigma);
    
    // Init the distance map
    Matf distanceMat;
    Transformation<double> T01(robotPose1.data());
    Transformation<double> T0Camera = T01.compose(TbaseCamera);
    EnhancedStereo stereo(Transformation<double>(),
                params.data(), params.data(), stereoParams);
    stereo.generatePlane(T0Camera.inverseCompose(TbasePlane), distanceMat,
         vector<Vector3d>{Vector3d(-0.1, -0.1, 0), Vector3d(-0.1 + 3 * 0.45, -0.1, 0),
                          Vector3d(-0.1 + 3 * 0.45, 0.5, 0), Vector3d(-0.1, 0.5, 0) } );
    imshow("distance", distanceMat);
    
    
    
    // Init the localizer
    PhotometricLocalization localizer(params.data(), params.data(), stereoParams);   
    localizer.initCloud(img1, distanceMat);
    
    
                    
    while (getline(paramFile, imageInfo))
    {
        istringstream imageStream(imageInfo);
        
        imageStream >> imageName;
        for (auto & x : robotPose2) imageStream >> x;
    
        Transformation<double> T02(robotPose2.data());
        Transformation<double> T12 = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
        cout << "REAL POSE : " << T12 << endl;
        T12 = T12.compose(Transformation<double>(-0.1, 0.05, -0.1, 0.03, 0.03, 0.05));
//        T12 = T12.compose(Transformation<double>(-0.001, 0.005, -0.01, 0.001, 0.001, 0.001));
        Matf img2 = imread(imageDir + imageName, 0);
        GaussianBlur(img2, img2, Size(0, 0), sigma, sigma);
        
        Matf img11, img12;
        localizer.wrapImage(img2, img11, T12);
        cout << T12 << endl;
        for (int iter = 0; iter < 1; iter++)
        {
            localizer.computePose(img2, T12);
            cout << T12 << endl;
            localizer.wrapImage(img2, img12, T12);
            imshow("img11", img11/256);
            imshow("img12", img12/256);
            imshow("delta img11", abs(img1 - img11)/250);
            imshow("delta img12", abs(img1 - img12)/250);
            waitKey(50);
        }
    }
    waitKey();
    return 0;
}
