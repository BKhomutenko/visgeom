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

void generateImages(const string baseName, const Transf xiBaseCam,
        const Transf xiOrigBoard, const vector<Transf> & odomVec)
{
    ofstream nameFile;
    nameFile.open(baseName + "_names.txt");
    
    array<double, 6> params {0.55, 1.1, 160, 160, 320, 240};
    EnhancedCamera camera(640, 480, params.data());
    
    BoardGenerator generator(&camera, 8, 5, 0.1);
    generator.setPlaneTransform(xiOrigBoard);    
//    Mat8u img0 = imread("/home/bogdan/projects/stack/aliasing4.png", 0);
    Mat8u img0 = imread("/home/bogdan/projects/stack/lena.png", 0);
    ImageGenerator generator2(&camera, img0, 100);
    generator2.setPlaneTransform(xiOrigBoard);
    
    Mat8u img;
    for (int i = 15; i < odomVec.size(); i++)
    {
        Transf xiCam = odomVec[i].compose(xiBaseCam);
        
        cout << xiCam << endl;
        
        generator2.generate(img, xiCam);
        
        const string name = baseName + to_string(i) + ".png";
        imwrite(name, img);
        nameFile << "            \"" + name + "\"," << endl;
    }
}

int main(int argc, char** argv) 
{
    ptree root;
    read_json(argv[1], root);
    vector<Transf> odom1Vec, odom2Vec;
    Transf xiBaseCam = readTransform(root.get_child("camera_transform"));
    Transf xiOrigBoard = readTransform(root.get_child("board_transform"));
    for (auto & odomItem : root.get_child("odom1"))
    {
        odom1Vec.emplace_back(readTransform(odomItem.second));
    }
    for (auto & odomItem : root.get_child("odom2"))
    {
        odom2Vec.emplace_back(readTransform(odomItem.second));
    }
    
    generateImages("lena_1_", xiBaseCam, xiOrigBoard, odom1Vec);
    generateImages("lena_2_", xiBaseCam, xiOrigBoard, odom2Vec);
    
    
    //stereo
    /*Transf xiBaseCam1(-0.35, 0, 0, 0, 0, 0);
    Transf xiBaseCam2(0.35, 0, 0, 0, 0, 0);
    Transf xiOrigBoard(-0.4, -0.2, 1.5, 0, 0, 0);
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            for (int k = -1; k < 2; k++)
            {
//                for (int l = -1; l < 2; l++)
//                {
                    const double step = 0.3;
                    odom1Vec.emplace_back(step * i, step * j, 0, 0, step * k, 0);
//                }
            }
        }
    }
    
    
    
    generateImages("images_stereo_1_", xiBaseCam1, xiOrigBoard, odom1Vec);
    generateImages("images_stereo_2_", xiBaseCam2, xiOrigBoard, odom1Vec);*/
    
    return 0;
}
