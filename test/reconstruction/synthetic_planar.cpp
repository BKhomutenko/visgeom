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
#include "projection/eucm.h"
#include "utils/image_generator.h"


int main(int argc, char** argv) 
{
    ptree root;
    read_json(argv[1], root);
    const vector<double> intrinsic = readVector(root.get_child("camera_intrinsics"));
    const int width = root.get<int>("image.width");
    const int heigth = root.get<int>("image.height");
    Transf xiCam0 = readTransform(root.get_child("camera_transform"));
    
    Mat8u foreImg = imread(root.get<string>("foreground"), 0);
    Mat8u backImg = imread(root.get<string>("background"), 0);
    
    EnhancedCamera camera(width, heigth, intrinsic.data());
    
    ImageGenerator generator(&camera, foreImg, 250);
//    generator.setBackground(backImg);
    const int iterMax = root.get<int>("steps");
    int cameraIncCount = 0;
    int boardPoseCount = 0;
    const string imageBaseName = root.get<string>("output_name");
    for (auto & boardPoseItem : root.get_child("plane_transform"))
    {
        generator.setPlaneTransform(readTransform(boardPoseItem.second));
        for (auto & cameraIncItem : root.get_child("camera_increment"))
        {
            Transf dxi = readTransform(cameraIncItem.second);
            Transf xiCam = xiCam0;
            
            for (int i = 0; i < iterMax; i++, xiCam = xiCam.compose(dxi))
            {
                cout << xiCam << endl;
                Mat8u dst;
                generator.generate(dst, xiCam);
                const string imgName = imageBaseName + "_" + to_string(boardPoseCount) + "_" +
                        to_string(cameraIncCount) + "_" + to_string(i) + ".png";
                cout << imgName << endl;
                imwrite(imgName, dst);
            }
            cameraIncCount++;
        }
        boardPoseCount++;
    }
    
    return 0;
}
