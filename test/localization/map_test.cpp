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

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    
    
    
    PhotometricMapping odom(root);
    
    
    string wordlFile = root.get<string>("render");
    ptree world;     
    read_json(wordlFile, world);
    RenderDevice device(world);
    device.setCamera(odom._camera);
    int depthMapCount = 0;
    for (auto & x : root.get_child("program"))
    {
        Transf xi = readTransform(x.second.get_child("initial"));
        odom.reInit(xi);
        for (auto & iteration : x.second.get_child("trajectory"))
        {
            Transf zeta = readTransform(iteration.second.get_child("increment"));
            int incrementCount = iteration.second.get<int>("step_count");
            
                for (int i = 0; i < incrementCount; i++, xi = xi.compose(zeta))
                {
                    device.setCameraTransform(xi.compose(odom._xiBaseCam));
                    Mat8u img0, img1;
            //        if (i > 0) odom.setDepth(device.getDepthBuffer());
//                    device.render(img0);
//                    cv::flip(img0, img1, 0);
                    
                    device.render(img1);
//                    Vector6d noise = Vector6d::Random() / 100;
//                    Transf eps(noise.data());
//                    odom.feedOdometry(xi.compose(eps));
                    odom.feedOdometry(xi);
                    odom.feedImage(img1);
                    cout << "GROUND TRUTH : " << endl;
                    cout << "    " << xi << endl;
                    cout << "ESTIMATE : " << endl;
                    cout << "    " << odom._interFrame.xi.compose(odom._xiLocal) << endl;
                    imshow("rendered", img1);
                    
                    Mat32f depthMat;
                    if (not odom._depth.empty())
                    { 
                        odom._depth.toMat(depthMat);
                        imshow("depth", depthMat / 10);
                        imwrite("depth" + to_string(depthMapCount++) + ".png", depthMat  * 25);
                    }
                    waitKey(50);
                }
    
        }
        
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

