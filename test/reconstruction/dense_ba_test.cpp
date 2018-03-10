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
#include "ceres.h"
#include "json.h"

#include "geometry/geometry.h"
#include "projection/eucm.h"
#include "projection/jacobian.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "render/render.h"
#include "localization/sparse_odom.h"
//FIXME make an argument
ofstream results;
string histDataName;


struct ReprojectionCost : ceres::SizedCostFunction<1, 6, 1>
{
    ReprojectionCost(const ICamera * camera,
            const Mat32f & img1, 
            const Mat32f & img2, 
            const int u, 
            const int v):
        _camera(camera->clone()), 
        _imageGrid(img2.cols, img2.rows, (float*)(img2.data)),
        _val(img1(v ,u))
    {
        Vector2d p(u, v);
        _camera->reconstructPoint(p, _dir);
        _dir.normalize();
    }

    virtual ~ReprojectionCost() { delete _camera; }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const
    {
        Transf xi(params[0]);
        Vector3d X1 = _dir * params[1][0];
        Vector3d X2;
        xi.inverseTransform(X1, X2);
        Vector2d pt;
        if (not _camera->projectPoint(X2, pt))
        {
            residual = 0;
            if (jacobian != NULL)
            {
                if (jacobian[0] != NULL)
                {
                    fill(jacobian[0], jacobian[0] + 6, 0);
                }
                if (jacobian[1] != NULL)
                {
                    jacobian[1][0] = 0;
                }
            }
            return false;
        }
        ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(_imageGrid);
        
        // image interpolation and gradient
        double f;
        Covector2d grad;
        imageInterpolator.Evaluate(pt[1], pt[0],
                &f, &grad[1], &grad[0]);
        
        *residual = f - _val;
        
        
        if (jacobian != NULL)
        {
            CameraJacobian jacobianCalculator(_camera, xi, Transf(0, 0, 0, 0, 0, 0));
            if (jacobian[0] != NULL)
            {
                Map<Covector6d> dfdxi(jacobian[0]);
                jacobianCalculator.dfdxi(X2, grad, dfdxi.data());
            }
            if (jacobian[1] != NULL)
            {
                Matrix23drm projJac;
                _camera->projectionJacobian(X2, projJac.data(), projJac.data() + 3);
                jacobian[1][0] = grad * projJac * xi.rotMatInv() * _dir;
            }
        }
        return true;
    }
    
    ICamera * _camera;
    Grid2D<float> _imageGrid;
    Vector3d _dir;
    float _val;
};


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
            ofs << (depthGT(vgt, ugt) - depth(v, u)) << "    " << depthGT(vgt, ugt) << "   " << sigma(v, u) << endl;
            if (sigma(v, u) > 10 or abs(depthGT(vgt, ugt) - depth(v, u)) >  sigma(v, u))
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
    
//    results << sqrt(err2 / N)*1000  << "    " << 100 * N / double(Nmax)
//        << "    " << dist / Nmax << "    " << N / double(Ngt);
     ofs.close();    
    imshow("inliers", inlierMat);
}

int main (int argc, char const* argv[])
{
    
    ptree root;
    read_json(argv[1], root);
    
    Transf xiCam = readTransform(root.get_child("trajectory.initial"));
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
    
    const string imageBaseName = "output_";
    
//    results.open("analysis_outoput.txt");
    
    int cameraIncCount = 0;
    //depth GT
    Mat32f depthGT, depth, sigmaMat;
    
    Mat8u img1;
    device.setCamera(&camera);
    device.setCameraTransform(xiCam);
    device.render(img1);
    depthGT = device.getDepthBuffer().clone();
    imshow("depthGT", depthGT / 10);
    
    
    
    
    
    
    
    
    
    // The initial approximation
    
    
    Transf err0 = Transf(0.00, -0.00, 0.00, 0.005, 0.004, 0.001);
    
    device.setCameraTransform(xiCam.compose(zeta));
    
    Mat8u img2;
    device.render(img2);
    
    
    SparseOdometry _sparseOdom(&camera, Transf(0, 0, 0, 0, 0, 0));
    
    //Compute the initial transformation
    _sparseOdom.feedData(img1, Transf(0, 0, 0, 0, 0, 0));
    _sparseOdom.feedData(img2, zeta.compose(err0));
    Transf xiStereo = _sparseOdom.getIncrement();      
    
    
//    cout << zeta.inverseCompose(xiStereo) << endl;
    
    EnhancedSgm stereo(zeta, &camera, &camera, stereoParams);
    
    DepthMap depthStereo;
    stereo.computeStereo(img1, img2, depthStereo);
    
    device.setCameraTransform(xiCam.compose(zeta).compose(zeta));
    Mat8u img3;
    device.render(img3);
    
    
    
    //Complete optimization
    
    array<double, 6> arr1, arr2;
    
    
    zeta.toArray(arr1.data());
    zeta.compose(zeta).compose(err0).toArray(arr2.data());
    
    Problem problem;
    
    Mat32f imgf1, imgf2, imgf3;
    img1.convertTo(imgf1, CV_32F);
    img2.convertTo(imgf2, CV_32F);
    img3.convertTo(imgf3, CV_32F);
    
    imgf1 /= 255;
    imgf2 /= 255;
    imgf3 /= 255;
    
    depthStereo.toMat(depth);
    depthStereo.sigmaToMat(sigmaMat);
    imshow("depth0", depth / 10);
    Mat_<double> depth64;
    
    depth.convertTo(depth64, CV_64F);
    
    for (int v = 0; v < img1.rows; v++)
    {
        for (int u = 0; u < img1.cols; u++)
        {
            if (depthStereo.at(u, v) == OUT_OF_RANGE) continue;
            problem.AddResidualBlock(new ReprojectionCost(&camera, imgf1, imgf2, u, v), NULL, arr1.data(), &depth64(v, u));
            problem.AddResidualBlock(new ReprojectionCost(&camera, imgf1, imgf3, u, v), NULL, arr2.data(), &depth64(v, u));
        }
    }
    
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 1000;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    cout << summary.FullReport() << endl;
    
    
    Transf zetaRes1(arr1.data());
    Transf zetaRes2(arr2.data());
    
    cout << err0 << endl;
    cout << zetaRes1.inverseCompose(zeta) << endl;
    cout << zetaRes1.inverseCompose(zetaRes2).inverseCompose(zeta) << endl;
    
    imshow("img1", img1);
    imshow("img2", img2);
    imshow("img3", img3);
    
    
    
    imshow("depth", depth64 / 10);
    
    
    analyzeError(depthGT, depth, sigmaMat, stereoParams);
    
    depth64.convertTo(depth, CV_64F);
    analyzeError(depthGT, depth, sigmaMat, stereoParams);
    
    waitKey();
    return 0;
}

int main_bak(int argc, char** argv) 
{

    ptree root;
    read_json(argv[1], root);
    
    Transf xiCam0 = readTransform(root.get_child("trajectory.initial"));
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
    
    const string imageBaseName = "output_";
    
    results.open("analysis_outoput.txt");
    
    int cameraIncCount = 0;
    //depth GT
    Mat32f depthGT, depth, sigmaMat;
    
    Mat8u img1;
    device.setCameraTransform(xiCam0);
    device.render(img1);
    depthGT = device.getDepthBuffer().clone();
    imshow("depthGT", depthGT / 10);
    //base frame
    const int noiseLevel = root.get<int>("noise");
    if (noiseLevel != 0)
    {
        Mat8u noise(img1.size());
        randu(noise, 0, noiseLevel);
        img1 -= noise;
    }
    //different increment diretion
    Transf xiCam = xiCam0.compose(zeta);
    cout << cameraIncCount << endl;
    
    MotionStereoParameters motionStereoParams(stereoParams);
    motionStereoParams.verbosity = 0;
    motionStereoParams.descLength = 7;
    
    MotionStereo motionStereo(&camera, &camera, motionStereoParams);
    motionStereo.setBaseImage(img1);
    //increment count
    DepthMap depthStereo;
    for (int i = 0; i < incrementCount; i++, xiCam = xiCam.compose(zeta))
    {
        histDataName = "hist_" + to_string(cameraIncCount) + "_" + to_string(i+1) + ".txt";
        std::ofstream ofs;
        ofs.open (histDataName, std::ofstream::out);
        ofs.close();
        
        device.setCameraTransform(xiCam);
        Mat8u img2;
        device.render(img2);
        if (noiseLevel != 0)
        {
            Mat8u noise(img2.size());
            randu(noise, 0, noiseLevel);
            img2 -= noise;
        }
        Transf TleftRight = xiCam0.inverseCompose(xiCam);
        cout << xiCam << endl;
        results << TleftRight.trans().norm() << "    ";
        
        if (i < 2 or root.get<bool>("sgm"))
        {
            EnhancedSgm stereo(TleftRight, &camera, &camera, stereoParams);
            stereo.computeStereo(img1, img2, depthStereo);
        }
        else
        {
            depthStereo = motionStereo.compute(TleftRight, img2, depthStereo);
        }
        
        
        if (root.get<bool>("filter_noise")) depthStereo.filterNoise();
        
        
        //output the point cloud
        if (i == 9)
        {
            MHPack pointPack;
            depthStereo.reconstruct(pointPack);
            ofstream cloudos("cloud.txt");
            
            Vector3dVec cloud;
            
            for (auto & x : pointPack.cloud)
            {
                cloudos << x.transpose() << endl;
            }
            cloudos.close();
        }
        
        
        
        
        depthStereo.toMat(depth);
        depthStereo.sigmaToMat(sigmaMat);
        analyzeError(depthGT, depth, sigmaMat, stereoParams);
        results << endl;
        imshow("depth", depth / 10);
        imshow("img", img2);
        imwrite("depth_" + to_string(i) + ".png", depth * 25);
        imwrite("img_" + to_string(i) + ".png", img2);
        waitKey();
    }
    results.close();
    return 0;
}



