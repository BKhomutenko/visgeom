// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"
#include "json.h"

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "reconstruction/depth_map.h"

#include "localization/mono_odom.h"

#include "render/render.h"

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    
    Transf xi = readTransform(root.get_child("trajectory.initial"));
    Transf zeta = readTransform(root.get_child("trajectory.increment"));
    int incrementCount = root.get<int>("trajectory.step_count");
    MonoOdometry odom(root);
    
    
    string wordlFile = root.get<string>("render");
    ptree world;     
    read_json(wordlFile, world);
    RenderDevice device(world);
    device.setCamera(odom._camera);
    
    for (int i = 0; i < incrementCount; i++, xi = xi.compose(zeta))
    {
        device.setCameraTransform(xi.compose(odom._xiBaseCam));
        Mat8u img1;
        device.render(img1);
        
        odom.feedWheelOdometry(xi);
        odom.feedImage(img1);
        
        Mat32f depthMat;
        odom.depth.toMat(depthMat);
        
        
        imshow("rendered", img1);
        imshow("depth", depthMat / 10);
        waitKey();
    }
    
    return 0;
}

