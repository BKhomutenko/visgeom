#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "timer.h"
#include "utils/curve_rasterizer.h"
#include "reconstruction/eucm_epipolar.h"
#include "reconstruction/epipolar_descriptor.h"
#include "reconstruction/depth_map.h"


Mat8u img1, img2;

int main(int argc, char** argv)
{	
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
    paramFile.ignore();
    Transformation<double> T01(robotPose1.data()), T02(robotPose2.data());
    
    Transformation<double> TleftRight = T01.compose(TbaseCamera).inverseCompose(T02.compose(TbaseCamera));
    
    EnhancedCamera cam1(params1.data()), cam2(params2.data());
    
    EnhancedEpipolar epipolar(&cam1, &cam2, TleftRight, 2000);
    
    string fileName1, fileName2;
    
    getline(paramFile, fileName1); //to SGM parameters
    
    const int LENGTH = 5;
    //TODO fix the constructor to avoid NULL 
    EpipolarDescriptor epipolarDescriptor(LENGTH, 3, {1, 2, 3, 5, 7, 9});
    StereoEpipoles epipoles(&cam1, &cam2, TleftRight);
    
    while(getline(paramFile, fileName1))
    {
        getline(paramFile, fileName2);
        
        img1 = imread(fileName1, 0);
        img2 = imread(fileName2, 0);
        
        Mat8u descStepMat(img1.size());
        
        for (int v = 0; v < img1.rows; v++)
        {
            for (int u = 0; u < img1.cols; u++)
            {
                Vector3d X;
                
                cam1.reconstructPoint(Vector2d(u, v), X);
                Vector2i pti(u, v);
                bool isInverted = epipoles.useInvertedEpipoleFirst(pti);
                CurveRasterizer<int, Polynomial2> raster(pti, epipoles.getFirstPx(isInverted),
                                         epipolar.getFirst(X));
                if (isInverted) raster.setStep(-1);
                vector<uint8_t> descriptor;
                const int step = epipolarDescriptor.compute(img1, raster, descriptor);
                if (step > 0)  descStepMat(v, u) = (10 - step) * 25; 
                else descStepMat(v, u) = 0;
            }
        }
        
        imshow("out1", img1);
        imshow("out2", img2);
        imshow("descStep", descStepMat);
        waitKey(); 
    }
    return 0;
}



