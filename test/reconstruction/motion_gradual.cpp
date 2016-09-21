#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"

#include "reconstruction/curve_rasterizer.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/eucm_sgm.h"

int main(int argc, char** argv)
{	
    /*Polynomial2 poly2;
    poly2.kuu = -1; 
    poly2.kuv = 1; 
    poly2.kvv= -1; 
    poly2.ku = 0.25; 
    poly2.kv = 0.25; 
    poly2.k1 = 5;
    
    CurveRasterizer<Polynomial2> raster(1, 1, -100, 100, poly2);
    CurveRasterizer2<Polynomial2> raster2(1, 1, -100, 100, poly2);

    auto tr0 = clock();
    int x1 = 0;
    int x2 = 0;
    for (int i = 0; i < 10000000; i++)
    {
        raster.step();
        x1 += raster.x;
    }
    auto tr1 = clock();
    
    for (int i = 0; i < 10000000; i++)
    {
        raster2.step();
        x2 += raster2.x;
    }
    auto tr2 = clock();
    
    cout << "optimized " << double(tr1 - tr0) / CLOCKS_PER_SEC << endl;
    cout << "simple " << double(tr2 - tr1) / CLOCKS_PER_SEC << endl;
    cout << x1 << " " << x2 << endl;
    return 0;*/
    ifstream paramFile(argv[1]);
    if (not paramFile.is_open())
    {
        cout << argv[1] << " : ERROR, file is not found" << endl;
        return 0;
    }
    
    array<double, 6> params;
    
    cout << "EU Camera model parameters :" << endl;
    for (auto & p: params) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 6> cameraPose;
    cout << "Camera pose wrt the robot :" << endl;
    for (auto & e: cameraPose) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbaseCamera(cameraPose.data());
    
    array<double, 6> robotPose1, robotPose2;

    
    
    MotionStereoParameters stereoParams;
    stereoParams.verbosity = 1;
    stereoParams.descLength = 5;
    
    SGMParameters stereoParams2;
    stereoParams2.salientPoints = false;
    stereoParams2.verbosity = 3;
    stereoParams2.hypMax = 1;
//    stereoParams.salientPoints = false;
    paramFile >> stereoParams2.u0;
    paramFile >> stereoParams2.v0;
    paramFile >> stereoParams2.dispMax;
    paramFile >> stereoParams2.scale;
    stereoParams.descLength = (stereoParams2.scale / 2 + 1) * 2 + 1;
    stereoParams.dispMax = 6;
    paramFile.ignore();
    string imageDir;
    getline(paramFile, imageDir);
    string imageInfo, imageName;
    getline(paramFile, imageInfo);
    istringstream imageStream(imageInfo);
    imageStream >> imageName;
    for (auto & x : robotPose1) imageStream >> x;

    Mat8u img1 = imread(imageDir + imageName, 0);
    cout << "Image name: "<< imageDir + imageName << endl;
    cout << "Image size: "<<img1.size()<<endl;;
    
    stereoParams2.u0 = 50;
    stereoParams2.v0 = 50;
    stereoParams2.uMax = img1.cols;
    stereoParams2.vMax = img1.rows;
    stereoParams2.setEqualMargin();
    
    int counter = 2;
    EnhancedCamera camera(params.data());
    DepthMap depth;
    depth.setDefault();
    
    MotionStereo stereo(&camera, &camera, stereoParams);
    stereo.setBaseImage(img1);
    
    
    //do SGM to init the depth
    getline(paramFile, imageInfo);
    getline(paramFile, imageInfo);
    getline(paramFile, imageInfo);
    imageStream.str(imageInfo);
    imageStream.clear();
    imageStream >> imageName;
    for (auto & x : robotPose2) imageStream >> x;
    Mat8u img2 = imread(imageDir + imageName, 0);
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    EnhancedSGM stereoSG(TleftRight, &camera, &camera, stereoParams2);
    
    Timer timer;
    stereoSG.computeStereo(img1, img2, depth);
    cout << timer.elapsed() << endl; 
    Mat32f res, sigmaRes;
    Mat32f res2, sigmaRes2;
    depth.toInverseMat(res2);
    depth.sigmaToMat(sigmaRes2);
    imshow("res" + to_string(counter), res2 *0.12);
    imshow("sigma " + to_string(counter), sigmaRes2*5);
    cv::waitKey(0);
    counter++;
    while (getline(paramFile, imageInfo))
    {
        istringstream imageStream(imageInfo);
        
        imageStream >> imageName;
        for (auto & x : robotPose2) imageStream >> x;
    
        Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
        Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
        
        Mat8u img2 = imread(imageDir + imageName, 0);

//        depth.setDefault();
        timer.reset();
        DepthMap depth2 = depth;
        stereo.computeDepth(TleftRight, img2, depth2);
        depth.merge(depth2);
        cout << timer.elapsed() << endl; 
        depth.toInverseMat(res);
        depth.sigmaToMat(sigmaRes);
//        imwrite(imageDir + "res" + to_string(counter++) + ".png", depth*200);
        imshow("res " + to_string(counter), res *0.12);
        imshow("sigma " + to_string(counter), sigmaRes*20);
        counter++; 
    }
    waitKey();
    return 0;
}



