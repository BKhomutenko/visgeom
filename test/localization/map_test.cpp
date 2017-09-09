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
    //FIXME temporary
//    Transf xiBaseCam( readTransform(root.get_child("xi_base_camera")) );
//    EnhancedCamera camera( readVector<double>(root.get_child("camera_params")).data() );
    
    string wordlFile = root.get<string>("render");
    ptree world;     
    read_json(wordlFile, world);
    RenderDevice device(world);
    device.setCamera(odom._camera); //FIXME
//    device.setCamera(&camera);
    int depthMapCount = 0;
    Mat8uc3 mapMat(Size(500, 500));
    mapMat.setTo(Scalar(255, 255, 255));
    int u0 = 250;
    int v0 = 50;
    int k = 30;
    ofstream fwo("wo.txt");
    ofstream fvo("vo.txt");
    ofstream fgt("gt.txt");
    int trajCount = 0;
    for (auto & x : root.get_child("program"))
    {
        Transf xi = readTransform(x.second.get_child("initial"));
        Transf xiwo = xi;
        
        //FIXME temporary
//        SparseOdometry odom(&camera, xiBaseCam);
    
        //FIXME temporary
        odom.reInit(xi);
        for (auto & iteration : x.second.get_child("trajectory"))
        {
            Transf zeta = readTransform(iteration.second.get_child("increment"));
            Transf zetawo = zeta;
            zetawo.rot()[1] += 0.005;
            int incrementCount = iteration.second.get<int>("step_count");
            
            for (int i = 0; i < incrementCount; i++, xi = xi.compose(zeta), xiwo = xiwo.compose(zetawo))
            {
                if (i % 2) continue;
                Vector6d noise = Vector6d::Random() * 0.002;
                Transf eps(noise.data());
                Transf xigt = xi.compose(eps);
//                device.setCameraTransform(xigt.compose(xiBaseCam));//FIXME
                device.setCameraTransform(xigt.compose(odom._xiBaseCam));
                Mat8u img0, img1;
            //        if (i > 0) odom.setDepth(device.getDepthBuffer());
//                    device.render(img0);
//                    cv::flip(img0, img1, 0);
                
                device.render(img1);
                odom.feedOdometry(xiwo); //FIXME
                
                
//                    if (not odom._depth.empty())
//                    { 
//                        Mat8u dst;
//                        odom._localizer.setDepth(odom._depth);
//                        odom._localizer.wrapImage(img1, dst, odom._xiLocal);
//                        Mat8u err;
//                        computeErr(odom._interFrame.img, dst, err);
//                        imshow("err_0", err);
//                    }
                
                cross(mapMat, xigt.trans()[0] * k + u0, xigt.trans()[2] * k + v0, 4, Scalar(255, 0, 0), 3);
                cross(mapMat, xiwo.trans()[0] * k + u0, xiwo.trans()[2] * k + v0, 3, Scalar(0, 255, 0));
                
                odom.feedImage(img1); ///FIXME
//                odom.feedData(img1, xiwo);
                

                Transf xivo =  odom._interFrame.xi.compose(odom._xiLocal); //FIXME
//                Transf xivo = odom.getIntegrated();




                cross(mapMat, xivo.trans()[0] * k + u0, xivo.trans()[2] * k + v0, 3, Scalar(0, 0, 255));
//                cout << "WHEEL ODOM : " << endl;
//                cout << "    " << xiwo << endl;
                cout << "GROUND TRUTH : " << endl;
                cout << "    " << xigt << endl;
                cout << "ESTIMATE : " << endl;
                cout << "    " << xivo << endl;
                fvo << xivo << endl;
                fwo << xiwo << endl;
                fgt << xigt << endl;
                //FIXME
//                cout << "ESTIMATE : " << endl;
//                cout << "    " << odom._interFrame.xi.compose(odom._xiLocal) << endl;
                imshow("rendered", img1);
                imshow("map", mapMat);
                
                
                
                
                //FIXME
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
        ofstream fkf("kf" + to_string(trajCount) + ".txt");
        for (auto & kf : odom._frameVec)
        {
            fkf << kf.xi << endl;
        }
        fkf.close();
        trajCount++;
    }
    fvo.close();
    fwo.close();
    fgt.close();
    
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

