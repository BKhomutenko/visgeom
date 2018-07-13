#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/depth_map.h"

EnhancedEpipolar * epipolar;
Mat8u img1, img2;
Mat8u orig1, orig2;
EnhancedCamera * cam1, * cam2;
EpipolarDescriptor * epipolarDescriptor;
StereoEpipoles * epipoles;
Transformation<double> TleftRight;
int epipolarLength;

void CallBackFunc1(int event, int x, int y, int flags, void* userdata)
{
    if  ( event == cv::EVENT_LBUTTONDOWN )
    {
        cout << "Left button of the mouse is clicked 1 - position (" << x << ", " << y << ")" << endl;
        
        Vector3d X;
        Vector2d pt;
        cam1->reconstructPoint(Vector2d(x, y), X);
        Vector2i pti(x, y);
        auto useInverted = epipoles->chooseEpipole(CAMERA_1, pti);
        CurveRasterizer<int, Polynomial2> raster(pti, epipoles->getPx(CAMERA_1, useInverted),
                             epipolar->get(CAMERA_1, X));
        if (useInverted) raster.setStep(-1);
        vector<uint8_t> descriptor;
        const int step = epipolarDescriptor->compute(orig1, raster, descriptor);
        cout << "step :" << step << endl;
        cout << "Descriptor :" << endl;
        for (auto & x : descriptor)
        {
            cout << setw(5) << int(x);
        }
        cout << endl;
        
        //get the sample sequence
        cam2->projectPoint(TleftRight.rotMatInv() * X, pt);
        Vector2i pti2 = round(pt);
        useInverted = epipoles->chooseEpipole(CAMERA_2, pti2);
        CurveRasterizer<int, Polynomial2> raster2(pti2, epipoles->getPx(CAMERA_2, useInverted),
                             epipolar->get(CAMERA_2, X));
        if (useInverted) raster2.setStep(-1);                     
        raster2.setStep(step);                     
        vector<uint8_t> samleVec;
        raster2.steps(-5);
        cout << "Samples :" << endl;
        for (int i = 0; i < 256; i++, raster2.step())
        {
            cout << setw(5) << int(orig2(raster2.v, raster2.u));
        }
        cout << endl;
        
        
        epipolar->traceEpipolarLine(x, y, img2, CAMERA_1, epipolarLength);
        cv::circle(img1, Point(x, y), 1, Scalar(128), -1);
        imshow("out1", img1);
        imshow("out2", img2);
    }
}

void CallBackFunc2(int event, int x, int y, int flags, void* userdata)
{
    if  ( event == cv::EVENT_LBUTTONDOWN )
    {
        cout << "Left button of the mouse is clicked 2 - position (" << x << ", " << y << ")" << endl;
        epipolar->traceEpipolarLine(x, y, img1,  CAMERA_2, epipolarLength);
        cv::circle(img2, Point(x, y), 1, Scalar(128), -1);
        imshow("out1", img1);
        imshow("out2", img2);
    }
}

int main(int argc, char** argv)
{	



    ptree root;
    read_json(argv[1], root);
    TleftRight = readTransform(root.get_child("stereo_transformation"));
    cam1 = new EnhancedCamera( readVector<double>(root.get_child("camera_params_left")).data() );
    cam2 = new EnhancedCamera( readVector<double>(root.get_child("camera_params_right")).data() );

    int length = root.get<int>("stereo_parameters.stereo_parameters.descriptor_size");
    int reps = root.get<int>("stereo_parameters.stereo_parameters.descriptor_response_thresh");
    epipolarLength = root.get<int>("stereo_parameters.stereo_parameters.disparity_max");;
    epipolarDescriptor = new EpipolarDescriptor(length, reps, {1});
    epipoles = new StereoEpipoles(cam1, cam2, TleftRight);
    epipolar = new EnhancedEpipolar(cam1, cam2, TleftRight, 2000);




   
    img1 = imread(root.get<string>("image_left"), 0);
    img2 = imread(root.get<string>("image_right"), 0);
    if(img1.empty()) cout << "Error in " << root.get<string>("image_left")<< endl;
    if(img2.empty()) cout << "Error in " << root.get<string>("image_right") << endl;
    img1.copyTo(orig1);
    img2.copyTo(orig2);
    imshow("out1", img1);
    imshow("out2", img2);
    cv::setMouseCallback("out1", CallBackFunc1, NULL);
    cv::setMouseCallback("out2", CallBackFunc2, NULL);
    waitKey();

    /*Polynomial2 poly2;
    poly2.kuu = -1; 
    poly2.kuv = 1; 
    poly2.kvv= -1; 
    poly2.ku = 0.25; 
    poly2.kv = 0.25; 
    poly2.k1 = 5;
    
    CurveRasterizer<Polynomial2> raster(1, 1, -100, 100, poly2);
    CurveRasterizer2<Polynomial2> raster2(1, 1, -100, 100, poly2);

    auto tr0 = clock();
    int x1 = 0;
    int x2 = 0;
    for (int i = 0; i < 10000000; i++)
    {
        raster.step();
        x1 += raster.x;
    }
    auto tr1 = clock();
    
    for (int i = 0; i < 10000000; i++)
    {
        raster2.step();
        x2 += raster2.x;
    }
    auto tr2 = clock();
    
    cout << "optimized " << double(tr1 - tr0) / CLOCKS_PER_SEC << endl;
    cout << "simple " << double(tr2 - tr1) / CLOCKS_PER_SEC << endl;
    cout << x1 << " " << x2 << endl;
    return 0;*/
    
    
    
    //Create two images and trace only epipolar lines
    
    
    
    
    
    /*
    //TODO NOT TO DELETE, MIGHT BE USEFUL IN FUTURE
    ifstream paramFile(argv[1]);
    if (not paramFile.is_open())
    {
        cout << argv[1] << " : ERROR, file is not found" << endl;
        return 0;
    }
    
    array<double, 6> params1;
    array<double, 6> params2;
    
    cout << "First EU Camera model parameters :" << endl;
    for (auto & p: params1) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    cout << "Second EU Camera model parameters :" << endl;
    for (auto & p: params2) 
    {
        paramFile >> p;
        cout << setw(10) << p;
    }
    cout << endl;
    paramFile.ignore();
    
    array<double, 6> cameraPose;
    cout << "Camera pose wrt the robot :" << endl;
    for (auto & e: cameraPose) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();
    Transformation<double> TbaseCamera(cameraPose.data());
    
    array<double, 6> robotPose1, robotPose2;
    cout << "First robot's pose :" << endl;
    for (auto & e: robotPose1) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    cout << "Second robot's pose :" << endl;
    for (auto & e: robotPose2) 
    {
        paramFile >> e;
        cout << setw(10) << e;
    }
    cout << endl;
    paramFile.ignore();paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    
    TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    const int LENGTH = 7;
    cam1 = new EnhancedCamera(params1.data());
    cam2 = new EnhancedCamera(params2.data());
    
    Vector2d pt;
    cout << cam1->projectPoint(TleftRight.trans(), pt) << endl;
    cout << pt << endl;
    cout << cam2->projectPoint(-TleftRight.transInv(), pt) << endl;
    cout << pt << endl;
    
    
    
    epipolarDescriptor = new EpipolarDescriptor(LENGTH, 10, {1, 2});
    epipoles = new StereoEpipoles(cam1, cam2, TleftRight);
    epipolar = new EnhancedEpipolar(cam1, cam2, TleftRight, 2000);
    
    string fileName1, fileName2;
    
    getline(paramFile, fileName1); //to SGM parameters
    cout << "Input:" << fileName1 << endl;
    
    while(getline(paramFile, fileName1))
    {
        cout << "Input:" << fileName1 << endl;
        getline(paramFile, fileName2);
        cout << "Input:" << fileName2 << endl;
        
        img1 = imread(fileName1, 0);
        img2 = imread(fileName2, 0);
        if(img1.empty()) cout << "Error in " << fileName1 << endl;
        if(img2.empty()) cout << "Error in " << fileName2 << endl;
        img1.copyTo(orig1);
        img2.copyTo(orig2);
        imshow("out1", img1);
        imshow("out2", img2);
        cv::setMouseCallback("out1", CallBackFunc1, NULL);
        cv::setMouseCallback("out2", CallBackFunc2, NULL);
        waitKey(); 
    }
    delete epipolar;
    
    */
    return 0;
}



