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

bool PhotometricLocalization::computeExtrinsic(const Matf & img1,
            const Matf & img2, const Matf & dist, Transformation<double> & xi)
{
    if (not initCloud(img1, dist)) return false;
    
    Problem problem;
    array<double, 6> pose = xi.toArray();
    typedef DynamicAutoDiffCostFunction<PhotometricError<EnhancedProjector>> photometricCF;
    PhotometricError<EnhancedProjector> * errorFunctor;
    errorFunctor = new PhotometricError<EnhancedProjector>(cam2.params, colorVec, cloud, img2);
    photometricCF * costFunction = new photometricCF(errorFunctor);
    costFunction->AddParameterBlock(6);
    costFunction->SetNumResiduals(cloud.size());
    problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), pose.data());
    
    //run the solver
    Solver::Options options;
    options.max_num_iterations = 500;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-10;
    options.parameter_tolerance = 1e-10;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    cout << summary.FullReport() << endl;
    
    xi = Transformation<double>(pose.data());
    return true;
}


bool PhotometricLocalization::initCloud(const cv::Mat_<float> & img1, const cv::Mat_<float> & dist)
{
    vector<Vector2d> imagePointVec;
    for (int v = 0; v < img1.rows; v++)
    {
        for (int u = 0; u < img1.cols; u++)
        {
            //TODO change it
            int ud = u, vd = v;
            
            if (dist(vd, ud) == 0) continue;
            colorVec.push_back(img1(v, u));
            imagePointVec.emplace_back(u, v);
        }
    }
    // TODO check the reconstruction and discard bad points
    cam1.reconstructPointCloud(imagePointVec, cloud);
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
