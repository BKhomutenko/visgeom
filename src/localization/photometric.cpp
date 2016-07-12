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

#include "localization/photometric.h"

#include "io.h"
#include "std.h"
#include "eigen.h"
#include "ocv.h"

#include <ceres/ceres.h>
#include <ceres/cubic_interpolation.h>

#include "camera/eucm.h"

using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::SoftLOneLoss;

bool PhotometricLocalization::computePose(const Mat32f & img2, Transformation<double> & T12)
{
    Problem problem;
    array<double, 6> pose = T12.toArray();
    typedef DynamicAutoDiffCostFunction<PhotometricError<EnhancedProjector>> photometricCF;
    PhotometricError<EnhancedProjector> * errorFunctor;
    errorFunctor = new PhotometricError<EnhancedProjector>(cam2.params, dataPack, img2);
    photometricCF * costFunction = new photometricCF(errorFunctor);
    costFunction->AddParameterBlock(6);
    costFunction->SetNumResiduals(dataPack.cloud.size());
    problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), pose.data());
    
    //run the solver
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 5;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    cout << summary.FullReport() << endl;
    
    T12 = Transformation<double>(pose.data());
    return true;
}

bool PhotometricLocalization::wrapImage(const Mat32f & img2, Mat32f & imgOut,
        Transformation<double> & T12)
{
    imgOut.create(cam1.height, cam1.width);
    imgOut.setTo(0);
    
    vector<Vector3d> cloud2;
    T12.inverseTransform(dataPack.cloud, cloud2);
    
    vector<Vector2d> pointVec2;
    cam2.projectPointCloud(cloud2, pointVec2);
    
    for (int i = 0; i < pointVec2.size(); i++)
    {
        int u = round(pointVec2[i][0]), v = round(pointVec2[i][1]);
        if (u < 0 or u >= img2.cols or v < 0 or v > img2.rows) continue;
        imgOut(dataPack.idxVec[i]) = img2(v, u);
    }
    return true;
}

const int step = 27;

bool PhotometricLocalization::initCloud(const Mat32f & img1, const Mat32f & dist)
{
    vector<Vector2d> imagePointVec;
    vector<double> distVec;
    for (int v = 0; v < img1.rows; v += step)
    {
        for (int u = 0; u < img1.cols; u += step)
        {
            //TODO change it
            int ud = params.uSmall(u);
            int vd = params.vSmall(v);
            if (ud < 0 or ud >= dist.cols or vd < 0 or vd >= dist.rows) continue;
            if (dist(vd, ud) < 0.01) continue;
//            cout << u << " " << v <<" " << ud << " " << vd << " " << dist(vd, ud) << endl;
            dataPack.colorVec.push_back(img1(v, u));
            imagePointVec.emplace_back(u, v);
            distVec.push_back(dist(vd, ud));
            dataPack.idxVec.push_back(v*img1.cols + u);
        }
    }
//    cout << indexVec.size() << endl;
    // TODO check the reconstruction and discard bad points
    cam1.reconstructPointCloud(imagePointVec, dataPack.cloud);
    for (int i = 0; i < dataPack.cloud.size(); i++)
    {
        dataPack.cloud[i] = dataPack.cloud[i] * (distVec[i] / dataPack.cloud[i].norm());
//        cout << cloud[i].transpose() << endl;
    }
    return true;
}

