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

#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "reconstruction/reproject_depth.h"
#include "reconstruction/eucm_stereo.h"

using namespace std;
using namespace cv;


int main(int argc, char** argv)
{	

//    ifstream paramFile(argv[1]);
//    if (not paramFile.is_open())
//    {
//        cout << argv[1] << " : ERROR, file is not found" << endl;
//        return 0;
//    }
//    
//    array<double, 6> params;
//    
//    cout << "EU Camera model parameters :" << endl;
//    for (auto & p: params) 
//    {
//        paramFile >> p;
//        cout << setw(15) << p;
//    }
//    cout << endl;
//    paramFile.ignore();
//    
//    array<double, 6> cameraPose;
//    cout << "Camera pose wrt the robot :" << endl;
//    for (auto & e: cameraPose) 
//    {
//        paramFile >> e;
//        cout << setw(15) << e;
//    }
//    cout << endl;
//    paramFile.ignore();
//    Transformation<double> TbaseCamera(cameraPose.data());
//    
//    array<double, 6> planePose;
//    cout << "Plane pose :" << endl;
//    for (auto & e: cameraPose) 
//    {
//        paramFile >> e;
//        cout << setw(15) << e;
//    }
//    cout << endl;
//    paramFile.ignore();
//    Transformation<double> TbasePlane(cameraPose.data());
//    
//    StereoParameters stereoParams;
//    paramFile >> stereoParams.uMargin;
//    paramFile >> stereoParams.vMargin;
//    paramFile >> stereoParams.dispMax;
//    paramFile >> stereoParams.scale;
//    
//    paramFile.ignore();
//    
//    string imageDir;
//    getline(paramFile, imageDir);
//    
//    string imageInfo, imageName;
//    array<double, 6> robotPose1, robotPose2;
//    getline(paramFile, imageInfo);
//    istringstream imageStream(imageInfo);
//    imageStream >> imageName;
//    for (auto & x : robotPose1) imageStream >> x;

//    Mat8u img1 = imread(imageDir + imageName, 0);
//    stereoParams.imageWidth = img1.cols;
//    stereoParams.imageHeight = img1.rows;

//    int counter = 2;
//    while (getline(paramFile, imageInfo))
//    {
//        istringstream imageStream(imageInfo);
//        
//        imageStream >> imageName;
//        for (auto & x : robotPose2) imageStream >> x;
//    
//        Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
//        Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
//        
//        cout << 111 << endl;
//        Mat8u img2 = imread(imageDir + imageName, 0);
//        cout << 111111 << endl;
//        EnhancedStereo stereo(TleftRight, params.data(), params.data(), stereoParams);

//        Mat8u res;
//        auto t2 = clock();
//         cout << 222 << endl;
//        stereo.comuteStereo(img1, img2, res);
//        auto t3 = clock();
////        cout << double(t3 - t2) / CLOCKS_PER_SEC << endl;
//        Mat32f distMat;
//        Mat32f planeMat;
//         cout << 333 << endl;
//        stereo.computeDistance(distMat);
//        Transformation<double> T0Camera = T01.compose(TbaseCamera);
//        stereo.generatePlane(T0Camera.inverseCompose(TbasePlane), planeMat,
//         Vector3dVec{Vector3d(-0.1, -0.1, 0), Vector3d(-0.1 + 3 * 0.45, -0.1, 0),
//                          Vector3d(-0.1 + 3 * 0.45, 0.5, 0), Vector3d(-0.1, 0.5, 0) } );
//        imshow("dist" + to_string(counter) , distMat/2);
//        imshow("disp" + to_string(counter) , res*4);
//        imshow("plane" , planeMat/2);
//        imwrite("/home/bogdan/projects/plane.png", planeMat);
//        double err = 0;
//        double err2 = 0;
//        double dist = 0;
//        int N = 0;
//        int Nmax = 0;
//        Mat32f inlierMat(planeMat.size());
//        inlierMat.setTo(0);
//        for (int u = 0; u < distMat.cols; u++)
//        {
//            for (int v = 0; v < distMat.rows; v++)
//            {
//                
//                if (planeMat(v, u) == 0) continue;
//                Nmax++;
//                dist += planeMat(v, u);
//                inlierMat(v, u) = 1;
//                if (distMat(v, u) == 0 or distMat(v, u) != distMat(v, u) or planeMat(v, u) != planeMat(v, u)) continue;
//                if (abs(distMat(v, u) - planeMat(v, u)) > 0.10) continue;
//                inlierMat(v, u) = 0;
//                err += distMat(v, u) - planeMat(v, u);
//                err2 += pow(distMat(v, u) - planeMat(v, u), 2);
//                N++;
//            }
//        }
////        cout << (counter - 1) * 7 << " & " << dist/ Nmax * 1000 << " & " << err / N *1000 << " & " << sqrt(err2 / N)*1000  
////                << " & " << 100 * N / double(Nmax) << "\\\\" << endl << "\\hline" << endl;
//        cout << "avg err : " << err / N *1000 << " avg err2 : " 
//<< sqrt(err2 / N)*1000  << " number of inliers : " << 100 * N / double(Nmax) << endl;
//        cout << planeMat.size() << " " <<  distMat.size() << endl;
//        imshow("diff" + to_string(counter), abs(planeMat - distMat));
//        imshow("inliers" + to_string(counter), inlierMat);
//        counter++;
//        
//    }
//    waitKey();
    return 0;
}



