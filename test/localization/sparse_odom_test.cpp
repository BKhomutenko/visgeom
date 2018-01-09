// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "localization/sparse_odom.h"
#include "projection/eucm.h"

#include "reconstruction/eucm_sgm.h"

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

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    Transf xiBaseCam = readTransform(root.get_child("xi_base_camera"));
    vector<Transf> odomVec;
    for (auto & x : root.get_child("odom"))
    {
        odomVec.emplace_back(readTransform(x.second));
    }
    vector<string> fileVec;
    for (auto & x : root.get_child("names"))
    {
        fileVec.emplace_back(x.second.get_value<string>());
    }
    
    assert(fileVec.size() == odomVec.size());
    
    
    vector<double> intrinsic = readVector<double>(root.get_child("camera_params"));
    EnhancedCamera camera(intrinsic.data());
//    
//    /************************/
//    cout << "Test the jacobian rank for 2 points" << endl;
//    
//    vector<Vector3d> xVec1;
//    
//    xVec1.emplace_back(1, 1, 5);
//    xVec1.emplace_back(-1, .5, 4);
//    
//    Transf xi12(0.1, 0.1, 0.1, 0.1, 0.1, 0.3);
//    
//    testJacobian(xVec1, xi12, &camera);
//    
//    /************************/
//    cout << "Test the jacobian rank for 3 points" << endl;
//    
//    xVec1.clear();
//    xVec1.emplace_back(1, 1, 5);
//    xVec1.emplace_back(-1, .5, 4);
//    xVec1.emplace_back(0, -1, 5);
//    
//    testJacobian(xVec1, xi12, &camera);
//    
//    /************************/
//    cout << "Test the jacobian rank for 4 points" << endl;
//    
//    xVec1.clear();
//    xVec1.emplace_back(1, 1, 5);
//    xVec1.emplace_back(-1, 1, 4);
//    xVec1.emplace_back(1, -1, 5);
//    xVec1.emplace_back(-1, -1, 4);
//    
//    testJacobian(xVec1, xi12, &camera);
//    
//    /************************/
//    cout << "Test the jacobian rank for 5 points" << endl;
//    
//    xVec1.clear();
//    xVec1.emplace_back(1, 1, 5);
//    xVec1.emplace_back(-1, 1, 4);
//    xVec1.emplace_back(1, -1, 5);
//    xVec1.emplace_back(-1, -1, 4);
//    xVec1.emplace_back(0, 0, 3);
//    
//    testJacobian(xVec1, xi12, &camera);
//    
//    /************************/
//    cout << "Test the jacobian rank for 6 points" << endl;
//    
//    xVec1.clear();
//    xVec1.emplace_back(1, 1, 5);
//    xVec1.emplace_back(-1, 1, 4);
//    xVec1.emplace_back(1, -1, 5);
//    xVec1.emplace_back(-1, -1, 4);
//    xVec1.emplace_back(-1, 0, 3);
//    xVec1.emplace_back(1, 0, 5);
//    
//    testJacobian(xVec1, xi12, &camera);
//    
    /****************************/
    cout << "odometry test" << endl;
    SparseOdometry odom(&camera, xiBaseCam);
    
    
    Transf xiBaseCam1(0.397979 , 1.55485 , 1.16112  ,  -1.59343, -0.00419144,   0.0167138);
    Transf xiCam12(0.00221874  , -0.300893, -0.00215483,  0.00839828,  0.00243455, -0.00758349);
    cout << xiBaseCam1.compose(xiCam12) << endl;
    //init stereoParameters
    SgmParameters stereoParams(root.get_child("stereo_parameters"));
    
    ofstream fodom("odom.txt");
    ofstream fsvo("svo_1.txt");
    Mat8u imgOld;
    for (int i = 0; i < fileVec.size(); i += 1)
    {
        cout << i << endl;
        Mat8u img = imread(fileVec[i], 0);
        odom.feedData(img, odomVec[i]);
        fodom << odomVec[0].inverseCompose(odomVec[i]) << endl;
        fsvo << odom.getIntegrated() << endl;
        Transf xi = xiBaseCam.inverseCompose(odom.getIncrement()).compose(xiBaseCam);
        
        
        if (not imgOld.empty())
        {
            EnhancedSgm sgm(xi, &camera, &camera, stereoParams);
            DepthMap depthStereo;
            sgm.computeStereo(imgOld, img, depthStereo);
            
            Mat32f depth;
            
            depthStereo.toMat(depth);
            imshow("depth", depth/55);
        }
        
        imshow("img", img);
        waitKey(50);
        img.copyTo(imgOld);
        
    }
    fodom.close();
    fsvo.close();
}
