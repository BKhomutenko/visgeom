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

#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
#include <iomanip>
#include <fstream>
#include <boost/format.hpp>

#include "geometry.h"
#include "vision.h"
#include "eucm.h"

using namespace cv;
using Eigen::Vector2d;
using Eigen::Vector3d;

typedef Mat_<float> fMat;

class Pinhole : public ICamera
{
public:

    Pinhole(double u0, double v0, double f)
    : ICamera(2*u0, 2*v0, 3) 
    {
        params[0] = u0;
        params[1] = v0;
        params[2] = f;
    }
    virtual ~Pinhole() {}

    virtual bool reconstructPoint(const Vector2d & src, Vector3d & dst) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & u = src(0);
        const double & v = src(1);
        dst << (u - u0)/f, (v - v0)/f, 1;
        return true;
    }

    /// projects 3D points onto the original image
    virtual bool projectPoint(const Vector3d & src, Vector2d & dst) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & x = src(0);
        const double & y = src(1);
        const double & z = src(2);
        if (z < 1e-2)
        {
            dst << -1, -1;
            return false;
        }
        dst << x * f / z + u0, y * f / z + v0;
        return true;
    }

    //TODO implement the projection and distortion Jacobian
    virtual bool projectionJacobian(const Vector3d & src, Eigen::Matrix<double, 2, 3> & Jac) const
    {
        const double & u0 = params[0];
        const double & v0 = params[1];
        const double & f = params[2];
        const double & x = src(0);
        const double & y = src(1);
        const double & z = src(2);
        double zz = z * z;
        Jac(0, 0) = f/z;
        Jac(0, 1) = 0;
        Jac(0, 2) = -x * f/ zz;
        Jac(1, 0)= 0;
        Jac(1, 1) = f/z;
        Jac(1, 2) = -y * f/ zz;
    }
    
    virtual Pinhole * clone() const
    {
        return new Pinhole(params[0], params[1], params[2]);
    }
};

void initRemap(const array<double, 6> & params1, const array<double, 3> & params2,
    fMat & mapX, fMat & mapY, const array<double, 3> & rot)
{
    EnhancedCamera cam1(params1.data());
    Pinhole cam2(params2[0], params2[1], params2[2]);
    vector<Vector2d> imagePoints;
    mapX.create(params2[1]*2, params2[0]*2);
    mapY.create(params2[1]*2, params2[0]*2);
    for (unsigned int i = 0; i < mapX.rows; i++)
    {
        for (unsigned int j = 0; j < mapX.cols; j++)
        {
            imagePoints.push_back(Vector2d(j, i));
        }
    }
    vector<Vector3d> pointCloud;
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
    
    fMat  mapX, mapY;
    initRemap(params, params2, mapX, mapY, rotation);
    
    string dirName, fileName;
    getline(paramFile, dirName);
    while (getline(paramFile, fileName))
    {
        fMat img = imread(dirName + fileName, 0);
        fMat img2;
        remap(img, img2, mapX, mapY, INTER_LINEAR);
        imshow("orig", img/255);
        imshow("res", img2/255);
        waitKey();
    }

    return 0;


}
