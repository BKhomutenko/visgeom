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
    paramFile >> stereoParams.u0;
    paramFile >> stereoParams.v0;
    paramFile >> stereoParams.disparityMax;
    paramFile >> stereoParams.blockSize;
    paramFile.ignore();
    
    string imageDir;
    getline(paramFile, imageDir);
    
    string imageInfo, imageName;
    array<double, 6> robotPose1, robotPose2;
    getline(paramFile, imageInfo);
    istringstream imageStream(imageInfo);
    imageStream >> imageName;
    for (auto & x : robotPose1) imageStream >> x;

    Matf img1 = imread(imageDir + imageName, 0);
    int counter = 2;
    while (getline(paramFile, imageInfo))
    {
        istringstream imageStream(imageInfo);
        
        imageStream >> imageName;
        for (auto & x : robotPose2) imageStream >> x;
    
        Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
        Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
        TleftRight = TleftRight.compose(Transformation<double>(0.01, 0.005, -0.001, 0.001, 0.001, 0.001));
        
        Matf img2 = imread(imageDir + imageName, 0);
        Matf planeMat;

        EnhancedStereo stereo(Transformation<double>(), img1.cols, img1.rows,
                params.data(), params.data(), stereoParams);

        Transformation<double> T0Camera = T01.compose(TbaseCamera);
        stereo.generatePlane(T0Camera.inverseCompose(TbasePlane), planeMat,
         vector<Vector3d>{Vector3d(-0.1, -0.1, 0), Vector3d(-0.1 + 3 * 0.45, -0.1, 0),
                          Vector3d(-0.1 + 3 * 0.45, 0.5, 0), Vector3d(-0.1, 0.5, 0) } );
       
        imshow("depth", planeMat);
        waitKey();
       
        PhotometricLocalization localizer(img1.cols, img1.rows,
                params.data(), params.data(), stereoParams);              
        
        cout << TleftRight << endl;
        localizer.computePose(img1, img2, planeMat, TleftRight);
        cout << TleftRight << endl;
    }
    return 0;
}
