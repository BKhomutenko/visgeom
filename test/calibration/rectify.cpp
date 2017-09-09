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
#include "json.h"

#include "geometry/geometry.h"
#include "projection/pinhole.h"
#include "projection/eucm.h"

void initRemap(const vector<double> & params1, const vector<double> & params2,
    Mat32f & mapX, Mat32f & mapY, const Transf & T)
{
    assert(params1.size() == 6);
    assert(params2.size() == 5);
    EnhancedCamera cam1(params1.data());
    Pinhole cam2(params2[2], params2[3], params2[4]);
    Vector2dVec imagePoints;
    mapX.create(params2[1], params2[0]);
    mapY.create(params2[1], params2[0]);
    for (unsigned int i = 0; i < mapX.rows; i++)
    {
        for (unsigned int j = 0; j < mapX.cols; j++)
        {
            imagePoints.push_back(Vector2d(j, i));
        }
    }
    Vector3dVec pointCloud;
    cam2.reconstructPointCloud(imagePoints, pointCloud);
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


    ptree root;
    read_json(argv[1], root);
    vector<double> paramsEucm = readVector<double>(root.get_child("camera_params"));

    //u0, v0, f
    vector<double> paramsPinhole = readVector<double>(root.get_child("pinhole_params"));
    Transf xi = readTransform(root.get_child("xi_eucm_pinhole"));
      
    Mat32f  mapX, mapY;
    initRemap(paramsEucm, paramsPinhole, mapX, mapY, xi);
    
    int count = 0;
    for (auto & x : root.get_child("image_names"))
    {
        Mat32f img = imread(x.second.get_value<string>(), 0);
        Mat32f img2;
        remap(img, img2, mapX, mapY, cv::INTER_LINEAR);
        imwrite("img_" + to_string(count++) + ".png", img2);
    }

    return 0;


}
