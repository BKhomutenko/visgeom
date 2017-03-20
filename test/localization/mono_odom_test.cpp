// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/depth_map.h"

#include "localization/photometric.h"
#include "localization/sparse_odom.h"

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    Transf xiBaseCam = readTransform(root.get_child("xiBaseCam"));
    vector<Transf> odomVec;
    for (auto & x : root.get_child("odometry"))
    {
        odomVec.emplace_back(readTransform(x.second));
    }
    vector<string> fileVec;
    for (auto & x : root.get_child("file_name"))
    {
        fileVec.emplace_back(x.second.get_value<string>());
    }
    
    assert(fileVec.size() == odomVec.size());
    
    string prefix = root.get<string>("prefix");
    
    vector<double> intrinsic = readVector(root.get_child("intrinsic"));
    EnhancedCamera camera(intrinsic.data());
    
    SparseOdometry odom(&camera, xiBaseCam);
    
    const int idx1 = 25, idx2 = 70, idx3 = 80;
    cout << prefix + fileVec[idx1] << endl;
    Mat8u img1 = imread(prefix + fileVec[idx1], 0);
    Mat8u img2 = imread(prefix + fileVec[idx2], 0);
    odom.feedData(img1, odomVec[idx1]);
    odom.feedData(img2, odomVec[idx2]);
    
    SGMParameters stereoParams;
    stereoParams.flawCost = 5;
    stereoParams.verbosity = 0;
    stereoParams.hypMax = 1;
//    stereoParams.salientPoints = false;
    stereoParams.u0 = 50;
    stereoParams.v0 = 50;
    stereoParams.dispMax = 26;
    stereoParams.scale = 2;
    
    stereoParams.uMax = img1.cols;
    stereoParams.vMax = img1.rows;
    stereoParams.setEqualMargin();
    cout << odom.getIncrement() << endl;
    Transf xiCam = xiBaseCam.inverseCompose(odom.getIncrement().compose(xiBaseCam));
    EnhancedSGM stereo(xiCam, &camera, &camera, stereoParams);
    
    DepthMap depth;
    Mat32f depthMat;
    stereo.computeStereo(img1, img2, depth);
    cout << "stereo's done" << endl;
    depth.toInverseMat(depthMat);
    imshow("out1", img2);
    imshow("res", depthMat);
    
    
    
    ScalePhotometric photometricLocalizer(5, &camera);
    photometricLocalizer.setVerbosity(1);
    
    photometricLocalizer.setXiBaseCam(xiBaseCam);
    imshow("img1", img1);
    waitKey();
    photometricLocalizer.computeBaseScaleSpace(img1);
    photometricLocalizer.setDepth(depth);
    
    MotionStereoParameters stereoMotionParams(stereoParams);
    
    MotionStereo mstereo(&camera, &camera, stereoMotionParams);
    mstereo.setBaseImage(img1);
    DepthMap d2 = depth;
    
    
    //TODO implement the whole localization loop
    for (int i = idx1; i < idx1 + 50; i++)
    {
        cout << "Idx : " << i << endl;
        Mat8u img3 = imread(prefix + fileVec[i], 0);
        Transf xi13(odomVec[idx1].inverseCompose(odomVec[i]));
        cout << xi13 << endl;
        photometricLocalizer.computePose(img3, xi13);
        cout << xi13 << endl; 
        if (i > idx1 + 20)
        {
            Transf xiCam12 = xiBaseCam.inverseCompose(xi13).compose(xiBaseCam);
            d2 = mstereo.compute(xiCam12, img3, d2);
//            d2.filterNoise();
            photometricLocalizer.setDepth(d2);
            d2.toInverseMat(depthMat); 
            //FIXME rewrite with Mat8u
            Mat32f img3f, img3fWrapped;
            img3.convertTo(img3f, CV_32F);
            photometricLocalizer.wrapImage(img3f, img3fWrapped, xi13);
            imshow("wrapped", img3fWrapped/256);
            imshow("res2", depthMat);
            
            //TODO for every new keyframe do SGBM
            
            DepthMap d3 = mstereo.compute(xiCam12, img3, d2, CAMERA_2);
//            d3.filterNoise();
            d3.toInverseMat(depthMat); 
            imshow("res3", depthMat);
            waitKey();
            waitKey(50);
        }
        
    }
    
    return 0;
}

