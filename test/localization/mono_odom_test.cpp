// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "reconstruction/eucm_sgm.h"
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
    
    const int idx1 = 0, idx2 = 25, idx3 = 25;
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
    depth.toInverseMat(depthMat);
    imshow("out1", img2);
    imshow("res", depthMat);
    
    Mat8u img3 = imread(prefix + fileVec[idx2], 0);
    
    ScalePhotometric photometricLocalizer(5, &camera);
    
    Transf xi23(xiBaseCam.inverseCompose(odomVec[idx2]).compose(xiBaseCam));
    photometricLocalizer.computeBaseScaleSpace(img1);
    photometricLocalizer.depth() = depth; //TODO make a method setDepth()
    photometricLocalizer.computePose(img3, xi23);
    cout <<xiBaseCam.compose(xi23).composeInverse(xiBaseCam) << endl;
    cout << odomVec[idx1].inverseCompose(odomVec[idx2]) << endl;
    waitKey();
}

