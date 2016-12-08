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

#include "geometry/geometry.h"
#include "projection/pinhole.h"
#include "projection/eucm.h"

void initRemap(const array<double, 6> & params1, const array<double, 3> & params2,
    Mat32f & mapX, Mat32f & mapY, const array<double, 3> & rot)
{
    EnhancedCamera cam1(params1.data());
    Pinhole cam2(params2[0], params2[1], params2[2]);
    Vector2dVec imagePoints;
    mapX.create(params2[1]*2, params2[0]*2);
    mapY.create(params2[1]*2, params2[0]*2);
    for (unsigned int i = 0; i < mapX.rows; i++)
    {
        for (unsigned int j = 0; j < mapX.cols; j++)
        {
            imagePoints.push_back(Vector2d(j, i));
        }
    }
    Vector3dVec pointCloud;
    cam2.reconstructPointCloud(imagePoints, pointCloud);
    Transformation<double> T(0, 0, 0, rot[0], rot[1], rot[2]);
    T.transform(pointCloud, pointCloud);
    cam1.projectPointCloud(pointCloud, imagePoints);
    
    auto pointIter = imagePoints.begin();
    for (unsigned int i = 0; i < mapX.rows; i++)
    {
        for (unsigned int j = 0; j < mapX.cols; j++)
        {
            mapX(i, j) = (*pointIter)[0];
            mapY(i, j) = (*pointIter)[1];
            ++pointIter;
        }
    }
}

int main(int argc, char** argv) {

    // read parameters and image names
    array<double, 6> params;
    ifstream paramFile(argv[1]);
    cout << argv[1] << endl;
    cout << "EU Camera model parameters :" << endl;
    for (auto & p: params) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 3> params2;
    cout << "Pinhole rectified camera parameters :" << endl;
    for (auto & p: params2) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 3> rotation;
    cout << "Rotation of the pinhole camera :" << endl;
    for (auto & p: rotation) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    Mat32f  mapX, mapY;
    initRemap(params, params2, mapX, mapY, rotation);
    
    string dirName, fileName;
    getline(paramFile, dirName);
    while (getline(paramFile, fileName))
    {
        Mat32f img = imread(dirName + fileName, 0);
        Mat32f img2;
        remap(img, img2, mapX, mapY, cv::INTER_LINEAR);
        imshow("orig", img/255);
        imshow("res", img2/255);
        waitKey();
    }

    return 0;


}
