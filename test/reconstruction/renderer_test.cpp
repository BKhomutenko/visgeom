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

string histDataName = "myHist.txt";
ofstream results;
void analyzeError(const Mat32f & depthGT, Mat32f & depth, 
        const Mat32f & sigma, const ScaleParameters & scaleParams)
{
    Mat8u inlierMat(depth.size());
    inlierMat.setTo(0);
    int Nmax = 0, Ngt = 0;
    double dist = 0;
    int N = 0;
    double err = 0, err2 = 0;
    double err3 = 0;
    std::ofstream ofs;
    ofs.open (histDataName, std::ofstream::out | std::ofstream::app);


 
    for (int u = 0; u < depth.cols; u++)
    {
        for (int v = 0; v < depth.rows; v++)
        {
            int ugt = scaleParams.uConv(u);
            int vgt = scaleParams.vConv(v);
            if (depthGT(vgt, ugt) != 0) Ngt++;
            
            if (depthGT(vgt, ugt) == 0 or depth(v, u) == 0 )
            {
                depth(v, u) = 0;
                continue;
                
            }
            
            if (depthGT(vgt, ugt) != depthGT(vgt, ugt) or depth(v, u) != depth(v, u)) continue;
            Nmax++;
            dist += depthGT(vgt, ugt);
            ofs << (depthGT(vgt, ugt) - depth(v, u)) / sigma(v, u) << endl;
            if (sigma(v, u) > 3 or abs(depthGT(vgt, ugt) - depth(v, u)) > 2.5 * sigma(v, u))
            {
                continue;
            }
            inlierMat(v, u) = 255;
            err += depthGT(vgt, ugt) - depth(v, u);
            err2 += pow(depthGT(vgt, ugt) - depth(v, u), 2);
            N++;
        }
    }
    cout << "avg err : " << err / N *1000 << " avg err2 : " 
        << sqrt(err2 / N)*1000  << " number of inliers : " << 100 * N / double(Nmax)
        << "   average distance : " << dist / Nmax << endl;
    
    results << sqrt(err2 / N)*1000  << "    " << 100 * N / double(Nmax)
        << "    " << dist / Nmax << "    " << N / double(Ngt);
     ofs.close();    
    imshow("inliers", inlierMat);
}

int main(int argc, char** argv) 
{

    
    ptree root;
    read_json(argv[1], root);
    vector<double> params = readVector<double>(root.get_child("camera_params"));
    RenderDevice device(root);
    EnhancedCamera camera(params.data());
    device.setCamera(&camera);
    
    Transf xi1(0., 0., 0., 0., 0., 0.);
    for (int i = 1; i < 10; i++)
    {
        
        
        Transf xi2(0.1 * i, 0., 0., 0., 0., 0.);
        
           
        Mat8u img1, img2;    
        device.setCameraTransform(xi1);
        device.render(img1);
        
        
//        
        device.setCameraTransform(xi2);
        device.render(img2);
        Mat32f depthGT = device.getDepthBuffer().clone();
        imshow("depth", device.getDepthBuffer() / 10);
         imshow("rendered", img2);
        waitKey();
    }
    
    
    return 0;
}
