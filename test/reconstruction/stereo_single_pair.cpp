#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "json.h"

#include "geometry/geometry.h"
#include "projection/eucm.h"
#include "reconstruction/eucm_sgm.h"
#include "reconstruction/depth_map.h"
#include "reconstruction/eucm_motion_stereo.h"

int main(int argc, char** argv)
{
    ptree root;
    read_json(argv[1], root);
    Transf zeta = readTransform(root.get_child("stereo_transformation"));
    EnhancedCamera camera_left( readVector<double>(root.get_child("camera_params_left")).data() );
    EnhancedCamera camera_right( readVector<double>(root.get_child("camera_params_right")).data() );

    //init stereoParameters
    SgmParameters stereoParams(root.get_child("stereo_parameters"));

    //the stereo is computed backwards
    EnhancedSgm sgm(zeta, &camera_left, &camera_right, stereoParams);

    //init depth and read images (change path and images names to correct)
    Mat32f depth;
    DepthMap depthStereo;
    Mat8u img1 = imread(root.get<string>("image_left"), 0);
    Mat8u img2 = imread(root.get<string>("image_right"), 0);

    //calc depth using SGM and convert it to Mat32f
    sgm.computeStereo(img1, img2, depthStereo);
    depthStereo.toInverseMat(depth);

    //visualisation (if necessary)
    imshow("depth", depth * (root.get<int>("brightness")/100.));
    imshow("img1", img1);
    imshow("img2", img2);
    waitKey();

}
