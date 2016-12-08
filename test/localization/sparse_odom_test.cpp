// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "utils/json.h"

#include "localization/sparse_odom.h"
#include "projection/eucm.h"

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
    
    for (int i = 0; i < fileVec.size(); i++)
    {
        cout << i << endl;
        Mat8u img = imread(prefix + fileVec[i], 0);
        odom.feedData(img, odomVec[i]);
        cout << odomVec[i] << endl;
        cout << odom.xiLocal << endl;
        imshow("img", img);
        
    }
    
}
