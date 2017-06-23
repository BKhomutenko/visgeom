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
#include "render/render.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"

//FIXME make an argument
ofstream results;

void analyzeError(const Mat32f & depthGT, const Mat32f & depth, 
        const Mat32f & sigma, const ScaleParameters & scaleParams)
{
    Mat8u inlierMat(depth.size());
    inlierMat.setTo(0);
    int Nmax = 0;
    double dist = 0;
    int N = 0;
    double err = 0, err2 = 0;
    for (int u = 0; u < depth.cols; u++)
    {
        for (int v = 0; v < depth.rows; v++)
        {
            int ugt = scaleParams.uConv(u);
            int vgt = scaleParams.vConv(v);
            if (depthGT(vgt, ugt) == 0 or depth(v, u) == 0 ) continue;
            Nmax++;
            dist += depthGT(vgt, ugt);
            if (depthGT(vgt, ugt) != depthGT(vgt, ugt) or depth(v, u) != depth(v, u)) continue;
            if (sigma(v, u) > 10 or abs(depthGT(vgt, ugt) - depth(v, u)) > 2.5 * sigma(v, u)) continue;
            inlierMat(v, u) = 255;
            err += depthGT(vgt, ugt) - depth(v, u);
            err2 += pow(depthGT(vgt, ugt) - depth(v, u), 2);
            N++;
        }
    }
    cout << "avg err : " << err / N *1000 << " avg err2 : " 
        << sqrt(err2 / N)*1000  << " number of inliers : " << 100 * N / double(Nmax)
        << "   average distance : " << dist / Nmax << endl;
    
//    results << sqrt(err2 / N)*1000  << "    " << 100 * N / double(Nmax)
//        << "    " << dist / Nmax;
        
    imshow("inliers", inlierMat);
}

int main(int argc, char** argv) 
{
    ptree root;
    read_json(argv[1], root);
    
    Transf xi = readTransform(root.get_child("trajectory.initial"));
    Transf zeta = readTransform(root.get_child("trajectory.increment"));
    int incrementCount = root.get<int>("trajectory.step_count");
    
    EnhancedCamera camera( readVector<double>(root.get_child("camera_params")).data() );
    
    
    
    //init stereoParameters
    SgmParameters stereoParams(root.get_child("stereo_parameters"));
    
    //the stereo is computed backwards
    EnhancedSgm sgm(zeta.inverse(), &camera, &camera, stereoParams);
    
    string wordlFile = root.get<string>("render");
    ptree world;     
    read_json(wordlFile, world);
    RenderDevice device(world);
    device.setCamera(&camera);
    
    Mat8u image0;
    Mat8u image1;
    
    for (int i = 0; i < incrementCount; i++, xi = xi.compose(zeta))
    {   
        device.setCameraTransform(xi);
        //depth GT
        Mat32f depthGT, depth, sigmaMat;
        
        image1.copyTo(image0);
        
        device.render(image1);
        
        if (i == 0) continue;
        
        device.getDepthBuffer().copyTo(depthGT);
        
        
        imshow("depthGT", depthGT / 10);
        
        DepthMap depthNew;
        sgm.computeStereo(image1, image0, depthNew);
        depthNew.toMat(depth);
        depthNew.sigmaToMat(sigmaMat);
        
        analyzeError(depthGT, depth, sigmaMat, depthNew);
        
        imshow("depthGT", depthGT / 10);
        imshow("depth", depth / 10);
        waitKey();
    }
    results.close();
    return 0;
}
