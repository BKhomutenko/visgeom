/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#include <vector>

#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>

#include <ceres/ceres.h>
#include <ceres/cubic_interpolation.h>

#include "localization/photometric.h"
#include "camera/eucm.h"

using std::vector;

using std::cout;
using std::endl;

using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::SoftLOneLoss;

typedef cv::Mat_<float> Matf;

bool PhotometricLocalization::computePose(const Matf & img2, Transformation<double> & T12)
{
    Problem problem;
    array<double, 6> pose = T12.toArray();
    typedef DynamicAutoDiffCostFunction<PhotometricError<EnhancedProjector>> photometricCF;
    PhotometricError<EnhancedProjector> * errorFunctor;
    errorFunctor = new PhotometricError<EnhancedProjector>(cam2.params, colorVec, cloud1, img2);
    photometricCF * costFunction = new photometricCF(errorFunctor);
    costFunction->AddParameterBlock(6);
    costFunction->SetNumResiduals(cloud1.size());
    problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), pose.data());
    
    //run the solver
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 300;
//    options.function_tolerance = 1e-10;
//    options.gradient_tolerance = 1e-10;
//    options.parameter_tolerance = 1e-10;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    cout << summary.FullReport() << endl;
    
    T12 = Transformation<double>(pose.data());
    return true;
}

bool PhotometricLocalization::wrapImage(const cv::Mat_<float> & img2, cv::Mat_<float> & imgOut,
        Transformation<double> & T12)
{
    imgOut.create(cam1.height, cam1.width);
    imgOut.setTo(0);
    
    vector<Vector3d> cloud2;
    T12.inverseTransform(cloud1, cloud2);
    
    vector<Vector2d> pointVec2;
    cam2.projectPointCloud(cloud2, pointVec2);
    
    for (int i = 0; i < pointVec2.size(); i++)
    {
        int u = round(pointVec2[i][0]), v = round(pointVec2[i][1]);
        if (u < 0 or u >= img2.cols or v < 0 or v > img2.rows) continue;
        imgOut(indexVec[i]) = img2(v, u);
    }
    return true;
}

const int step = 15;

bool PhotometricLocalization::initCloud(const cv::Mat_<float> & img1, const cv::Mat_<float> & dist)
{
    vector<Vector2d> imagePointVec;
    vector<double> distVec;
    for (int v = 0; v < img1.rows; v += step)
    {
        for (int u = 0; u < img1.cols; u += step)
        {
            //TODO change it
            int ud = (u - u0) / blockSize;
            int vd = (v - v0) / blockSize;
            if (ud < 0 or ud >= dist.cols or vd < 0 or vd >= dist.rows) continue;
            if (dist(vd, ud) < 0.01) continue;
//            cout << u << " " << v <<" " << ud << " " << vd << " " << dist(vd, ud) << endl;
            colorVec.push_back(img1(v, u));
            imagePointVec.emplace_back(u, v);
            distVec.push_back(dist(vd, ud));
            indexVec.push_back(v*img1.cols + u);
        }
    }
//    cout << indexVec.size() << endl;
    // TODO check the reconstruction and discard bad points
    cam1.reconstructPointCloud(imagePointVec, cloud1);
    for (int i = 0; i < cloud1.size(); i++)
    {
        cloud1[i] = cloud1[i] * (distVec[i] / cloud1[i].norm());
//        cout << cloud[i].transpose() << endl;
    }
    return true;
}
/*
class PhotometricLocalization
{
public:
    
    Transformation<double> computeExtrinsic(const cv::Mat_<float> & img1,
            const cv::Mat_<float> & img2,
            const cv::Mat_<float> & dist);
    
    void initCloud(const cv::Mat_<float> & dist);
    
private:
    EnhancedCamera cam1, cam2;
    vector<Vector3d> cloud;
    vector<float> colorVec;
};
*/
