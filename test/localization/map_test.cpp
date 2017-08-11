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
                    if (i > 2)
                    {
                        Vector6d noise = Vector6d::Random() / 100;
                        Transf eps(noise.data());
                        odom.feedOdometry(xi.compose(eps));
                    }
                    else
                    {
                        odom.feedOdometry(xi);
                    }
                    
//                    if (not odom._depth.empty())
//                    { 
//                        Mat8u dst;
//                        odom._localizer.setDepth(odom._depth);
//                        odom._localizer.wrapImage(img1, dst, odom._xiLocal);
//                        Mat8u err;
//                        computeErr(odom._interFrame.img, dst, err);
//                        imshow("err_0", err);
//                    }
                    
                    odom.feedImage(img1);
                    cout << "GROUND TRUTH : " << endl;
                    cout << "    " << xi << endl;
                    cout << "ESTIMATE : " << endl;
                    cout << "    " << odom._interFrame.xi.compose(odom._xiLocal) << endl;
                    imshow("rendered", img1);
                    
                    
                    
                    
                    
                    
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
                        imshow("depth", depthMat / 10);
//                        imwrite("depth" + to_string(depthMapCount++) + ".png", depthMat  * 25);
                    }
                    if (xi.trans().norm() > 20.65) waitKey();
                    else waitKey(50);
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

