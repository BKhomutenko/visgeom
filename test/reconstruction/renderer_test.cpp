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
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

#include "render/render.h"
#include "render/background.h"
#include "render/plane.h"

int main(int argc, char** argv) 
{


    ptree root;
    read_json(argv[1], root);
    array<double, 6> params {0.5, 1, 150, 150, 250, 250};
    Renderer renderer(root);
    
    renderer.setCameraParams(params.data());
    for (int i = 0; i < 100; i++)
    {
        renderer.setCameraTransform(Transf(0, 0, i*0.007, i*0.001, 0, -i*0.001));
//        renderer.setCameraTransform(Transf());
        cout << "initialized" << endl;
        renderer.fillBuffers();
        cout << "buffers" << endl;
        Mat8u res;
        renderer.fillImage(res);
        cout << "image" << endl;
        imshow("depth", renderer._depthMat/5);
        imshow("res", res);
        waitKey(1);
    }
    waitKey();
    return 0;
}
