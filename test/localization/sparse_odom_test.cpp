// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

//#include "localization/sparse_odom.h"
//#include "projection/eucm.h"

//#include "reconstruction/eucm_sgm.h"

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/depth_map.h"

#include "localization/mapping.h"

//To test the jacobian rank
#include "localization/local_cost_functions.h"
#include "reconstruction/triangulator.h"

//void testJacobian(const vector<Vector3d> & xVec1, Transf xi12,
//        const ICamera * camera)
//{
//    
//    vector<Vector2d> pVec2;
//    vector<Vector3d> xVec2;
//    
//    xi12.inverseTransform(xVec1, xVec2);
//    camera->projectPointCloud(xVec2, pVec2);
//    
//    SparseReprojectCost * projectionCost = new SparseReprojectCost(camera,
//                                                    xVec1, xVec2, pVec2, Transf());
//    
//    Triangulator triangulator(xi12);
//    vector<double> lambdaVec(xVec1.size());
//    triangulator.computeRegular(xVec1, xVec2, lambdaVec.data(), NULL, NULL, NULL);
//    for (auto & x : lambdaVec)
//    {
//        cout << x << endl;
//    }
//    array<double, 6> transfArr(xi12.toArray());
//    Matrix<double, Eigen::Dynamic, 1> residual(xVec1.size() * 2);
//    using MyJacMat = Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//    MyJacMat jacMat((xVec1.size() * 2), 6);
//    double * param = transfArr.data(), * jac = jacMat.data();
//    projectionCost->Evaluate(&param, residual.data(), &jac);
//    JacobiSVD<MyJacMat> svd(jacMat, Eigen::ComputeThinV);
//    cout << xi12 << endl;
//    cout << svd.singularValues() << endl;
//    cout << svd.matrixV() << endl;
//    cout << residual << endl;
//}

//int main(int argc, char** argv)
//{
//    ptree root;
//    read_json(argv[1], root);
//    Transf xiBaseCam = readTransform(root.get_child("xi_base_camera"));
//    vector<Transf> odomVec;
//    for (auto & x : root.get_child("odom"))
//    {
//        odomVec.emplace_back(readTransform(x.second));
//    }
//    vector<string> fileVec;
//    for (auto & x : root.get_child("names"))
//    {
//        fileVec.emplace_back(x.second.get_value<string>());
//    }
//    
//    assert(fileVec.size() == odomVec.size());
//    
//    
//    vector<double> intrinsic = readVector<double>(root.get_child("camera_params"));
//    EnhancedCamera camera(intrinsic.data());
////    
////    /************************/
////    cout << "Test the jacobian rank for 2 points" << endl;
////    
////    vector<Vector3d> xVec1;
////    
////    xVec1.emplace_back(1, 1, 5);
////    xVec1.emplace_back(-1, .5, 4);
////    
////    Transf xi12(0.1, 0.1, 0.1, 0.1, 0.1, 0.3);
////    
////    testJacobian(xVec1, xi12, &camera);
////    
////    /************************/
////    cout << "Test the jacobian rank for 3 points" << endl;
////    
////    xVec1.clear();
////    xVec1.emplace_back(1, 1, 5);
////    xVec1.emplace_back(-1, .5, 4);
////    xVec1.emplace_back(0, -1, 5);
////    
////    testJacobian(xVec1, xi12, &camera);
////    
////    /************************/
////    cout << "Test the jacobian rank for 4 points" << endl;
////    
////    xVec1.clear();
////    xVec1.emplace_back(1, 1, 5);
////    xVec1.emplace_back(-1, 1, 4);
////    xVec1.emplace_back(1, -1, 5);
////    xVec1.emplace_back(-1, -1, 4);
////    
////    testJacobian(xVec1, xi12, &camera);
////    
////    /************************/
////    cout << "Test the jacobian rank for 5 points" << endl;
////    
////    xVec1.clear();
////    xVec1.emplace_back(1, 1, 5);
////    xVec1.emplace_back(-1, 1, 4);
////    xVec1.emplace_back(1, -1, 5);
////    xVec1.emplace_back(-1, -1, 4);
////    xVec1.emplace_back(0, 0, 3);
////    
////    testJacobian(xVec1, xi12, &camera);
////    
////    /************************/
////    cout << "Test the jacobian rank for 6 points" << endl;
////    
////    xVec1.clear();
////    xVec1.emplace_back(1, 1, 5);
////    xVec1.emplace_back(-1, 1, 4);
////    xVec1.emplace_back(1, -1, 5);
////    xVec1.emplace_back(-1, -1, 4);
////    xVec1.emplace_back(-1, 0, 3);
////    xVec1.emplace_back(1, 0, 5);
////    
////    testJacobian(xVec1, xi12, &camera);
////    
//    /****************************/
//    cout << "odometry test" << endl;
//    SparseOdometry odom(&camera, xiBaseCam);
//    
//    
//    Transf xiBaseCam1(0.397979 , 1.55485 , 1.16112  ,  -1.59343, -0.00419144,   0.0167138);
//    Transf xiCam12(0.00221874  , -0.300893, -0.00215483,  0.00839828,  0.00243455, -0.00758349);
//    cout << xiBaseCam1.compose(xiCam12) << endl;
//    //init stereoParameters
//    SgmParameters stereoParams(root.get_child("stereo_parameters"));
//    
//    ofstream fodom("odom.txt");
//    ofstream fsvo("svo_1.txt");
//    Mat8u imgOld;
//    for (int i = 0; i < fileVec.size(); i += 1)
//    {
//        cout << i << endl;
//        Mat8u img = imread(fileVec[i], 0);
//        odom.feedData(img, odomVec[i]);
//        fodom << odomVec[0].inverseCompose(odomVec[i]) << endl;
//        fsvo << odom.getIntegrated() << endl;
//        Transf xi = xiBaseCam.inverseCompose(odom.getIncrement()).compose(xiBaseCam);
//        
//        
//        if (not imgOld.empty())
//        {
//            EnhancedSgm sgm(xi, &camera, &camera, stereoParams);
//            DepthMap depthStereo;
//            sgm.computeStereo(imgOld, img, depthStereo);
//            
//            Mat32f depth;
//            
//            depthStereo.toMat(depth);
//            imshow("depth", depth/55);
//        }
//        
//        imshow("img", img);
//        waitKey(50);
//        img.copyTo(imgOld);
//        
//    }
//    fodom.close();
//    fsvo.close();
//}

int main(int argc, char** argv)
{
    
    ptree root;
    read_json(argv[1], root);
    PhotometricMapping odom(root);
   
    Transf xiMap;
    bool xiMapInit = false;
    // construct the map
    int mapCount = 0;
    for (auto & fileName : root.get_child("map_files"))
    {
        mapCount++;
        ptree trajectoryData;
        read_json(fileName.second.get_value<string>(), trajectoryData);
        
        vector<Transf> gtVec;
        for (auto & x : trajectoryData.get_child("gt"))
        {
            gtVec.emplace_back(readTransform(x.second));
        }
        
        vector<string> fnameVec;
        for (auto & x : trajectoryData.get_child("names"))
        {
            fnameVec.emplace_back(x.second.get_value<string>());
        }
        if (not xiMapInit)
        {
            xiMapInit = true;
            xiMap = gtVec[0];
        }
        ofstream fgt("gtMap" + to_string(mapCount) + ".txt");
        for (int i = 0; i < fnameVec.size(); i++)
        {
            Mat8u img1 = imread(fnameVec[i], 0);
            if(odom.constructMap(xiMap.inverseCompose(gtVec[i]), img1))
            {
                fgt << xiMap.inverseCompose(gtVec[i]) << endl;
            }
        }
        fgt.close();
        continue;
    }
    
    // localize    
    int trajCount = 0;
    for (auto & fileName : root.get_child("trajectory_files"))
    {
        trajCount++;
        
        ptree trajectoryData;
        read_json(fileName.second.get_value<string>(), trajectoryData);
        
        vector<Transf> odomVec;
        for (auto & x : trajectoryData.get_child("odom"))
        {
            odomVec.emplace_back(readTransform(x.second));
        }
        
        vector<Transf> gtVec;
        for (auto & x : trajectoryData.get_child("gt"))
        {
            gtVec.emplace_back(readTransform(x.second));
        }
        
        vector<string> fnameVec;
        for (auto & x : trajectoryData.get_child("names"))
        {
            fnameVec.emplace_back(x.second.get_value<string>());
        }
        
        ofstream fvo("vo" + to_string(trajCount) + ".txt");
        ofstream fwo("wo" + to_string(trajCount) + ".txt");
        ofstream fkf("kf" + to_string(trajCount) + ".txt");
        ofstream fgt("gt" + to_string(trajCount) + ".txt");
        
        
        bool initialized = false;
        for (int i = 0; i < fnameVec.size(); i++)
        {
            //find the starting point
            if (not initialized)
            {
                Transf xi = xiMap.inverseCompose(gtVec[i]); //initial position
                if (odom.selectMapFrame(xi) == -1) continue;                
                cout << xiMap << endl << gtVec[i] << endl;
                odom.reInit(xi);
                initialized = true;
            }
            
            Mat8u img1 = imread(fnameVec[i], 0);
            odom._sparseOdom.feedData(img1, odomVec[i]);
            
            fvo << odom._sparseOdom.getIntegrated() << endl;
//            fwo << xiMap.inverseCompose(odomVec[i]) << endl; //true only if odomVec contains GT
            fwo << odomVec[i] << endl;
            fgt << xiMap.inverseCompose(gtVec[i]) << endl;

        }
        fvo.close();
        fwo.close();
        fkf.close();
        fgt.close();
    }
    
    
//    for (int i = 0; i < incrementCount; i++, xi = xi.compose(zeta))
//    {
//        device.setCameraTransform(xi.compose(odom._xiBaseCam));
//        Mat8u img1;
////        if (i > 0) odom.setDepth(device.getDepthBuffer());
//        device.render(img1);
//        
//        
//        Vector6d noise = Vector6d::Random() / 100;
//        Transf eps(noise.data());
//        odom.feedOdometry(xi.compose(eps));
////        odom.feedOdometry(xi);
//        odom.feedImage(img1);
//        cout << "GROUND TRUTH : " << endl;
//        cout << "    " << xi << endl;
//        cout << "ESTIMATE : " << endl;
//        cout << "    " << odom._interFrame.xi.compose(odom._xiLocal) << endl;
//        imshow("rendered", img1);
//        
//        Mat32f depthMat;
//        if (not odom._depth.empty())
//        { 
//            odom._depth.toMat(depthMat);
//            imshow("depth", depthMat / 10);
//        }
//        waitKey();
//    }
    
    return 0;
}
