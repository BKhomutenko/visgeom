// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "localization/mono_odom.h"

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
    MotionStereoParameters params;
    
    //TODO make reconstruction/utils and put there parameters readout
    params.scale = root.get<int>("scale");
    params.dispMax = root.get<int>("disparity_max");
    params.u0 = root.get<int>("u0");
    params.v0 = root.get<int>("v0");
    params.uMax = root.get<int>("u_max");
    params.vMax = root.get<int>("v_max");
    if (root.get<bool>("equal margins")) params.setEqualMargin();
    else
    {
        params.setXMargin(root.get<int>("u_marg"));
        params.setYMargin(root.get<int>("v_marg"));
    }
    
    MonoOdometry odom(&camera, xiBaseCam, params);
    
    for (int i = 0; i < fileVec.size(); i++)
    {
        cout << i << endl;
        Mat8u img = imread(prefix + fileVec[i], 0);
        odom.feedData(img, odomVec[i]);
        Mat32f depth;
        odom.depthMap.toInverseMat(depth);
        cout << odomVec[i] << endl;
        cout << odom._xiLocal << endl;
        imshow("img", img);
        imshow("depth", depth / 10);
        waitKey();
    }
    
}
