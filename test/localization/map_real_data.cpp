// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/depth_map.h"

#include "localization/mapping.h"

#include "render/render.h"

void computeErr(const Mat8u & img1, const Mat8u & img2, Mat8u & dst)
{
    dst.create(img1.size());
    for (int v = 0; v < img1.rows; v++)
    {
        for (int u = 0; u < img1.cols; u++)
        {
            if (img1(v, u) == 0 or img2(v, u) == 0)
            {
                dst(v, u) = 0;
                continue;
            }
            dst(v, u) = 127 + img1(v, u) - img2(v, u);
        }
    }
}

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    PhotometricMapping odom(root);
    int depthMapCount = 0;
    int trajCount = 0;
    for (auto & fileName : root.get_child("trajectory_files"))
    {
        trajCount++;
        ptree trajectoryData;
        read_json(fileName.second.get_value<string>(), trajectoryData);
        
        ofstream fvo("vo" + to_string(trajCount) + ".txt");
        ofstream fwo("wo" + to_string(trajCount) + ".txt");
        
        vector<Transf> odomVec;
        for (auto & x : trajectoryData.get_child("odom"))
        {
            odomVec.emplace_back(readTransform(x.second));
        }
        
        vector<string> fnameVec;
        for (auto & x : trajectoryData.get_child("names"))
        {
            fnameVec.emplace_back(x.second.get_value<string>());
        }
        
        Transf xi(0, 0, 0, 0, 0, 0); //initial position
        odom.reInit(xi);
        for (int i = 0; i < fnameVec.size(); i++)
        {
            odom.feedOdometry(odomVec[i]);
            Mat8u img1 = imread(fnameVec[i], 0);
            odom.feedImage(img1);
            
            
            
            cout << "ODOMETRY : " << endl;
            cout << "    " << odomVec[0].inverseCompose(odomVec[i]) << endl;
            cout << "ESTIMATE : " << endl;
            cout << "    " << odom._interFrame.xi.compose(odom._xiLocal) << endl;
            fvo << odom._interFrame.xi.compose(odom._xiLocal) << endl;
            fwo << odomVec[0].inverseCompose(odomVec[i]) << endl;
            if (not odom._depth.empty())
            {   
                Mat8u dst;
                Mat32f depthMat;
                odom._localizer.setDepth(odom._depth);
                odom._localizer.wrapImage(img1, dst, odom._xiLocal);
                imshow("wrapped", dst);
                imshow("base", odom._interFrame.img);
                Mat8u err;
                computeErr(odom._interFrame.img, dst, err);
                imshow("err", err);
                
                odom._depth.toMat(depthMat);
                imshow("depth", depthMat / 100);
                waitKey(50);
//                        imwrite("depth" + to_string(depthMapCount++) + ".png", depthMat  * 25);
            }

        }
        fvo.close();
        fwo.close();
        ofstream fkf("kf" + to_string(trajCount) + ".txt");
        for (auto & kf : odom._frameVec)
        {
            fkf << kf.xi << endl;
        }
        fkf.close();
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

