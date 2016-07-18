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

void ScalePhotometric::computeBaseScaleSpace(const Mat32f & img1)
{
    if (verbosity > 0) cout << "ScalePhotometric::computeBaseScaleSpace" << endl;
    scaleSpace1.generate(img1);
}

PhotometricPack ScalePhotometric::initPhotometricData(int scaleIdx)
{
    if (verbosity > 2) cout << "ScalePhotometric::initPhotometricData" << endl;
    PhotometricPack dataPack;
    dataPack.scaleIdx = scaleIdx;
    scaleSpace1.setActiveScale(scaleIdx);
    const Mat32f & img1 = scaleSpace1.get();
    const Mat32f & gradU1 = scaleSpace1.getGradU();
    const Mat32f & gradV1 = scaleSpace1.getGradV();
    double scale = scaleSpace1.getActiveScale();
    vector<double> distVec;
    vector<Vector2d> imagePointVec;
    vector<Vector2d> imageGradPointVec;
    if (verbosity > 3) cout << "    scaled image size : " << img1.size() << endl;
    for (int vs = 0; vs < img1.rows; vs++)
    {
        for (int us = 0; us < img1.cols; us++)
        {
            double gu = gradU1(vs, us);
            double gv = gradV1(vs, us);
            if (gu*gu + gv*gv < GRAD_THRESH) continue;
            int ub = scaleSpace1.uBase(us);
            int vb = scaleSpace1.vBase(vs);
            double dist = depthMap.nearest(ub, vb);
            if (dist < 0.01) continue;
            if (verbosity > 3) cout << "    " << vs << " " << us << endl;
            dataPack.colorVec.push_back(img1(vs, us));
            dataPack.gradValVec.push_back(sqrt(gu*gu + gv*gv) / scale);
            imagePointVec.emplace_back(ub, vb);
            imageGradPointVec.emplace_back(ub + gu/GRAD_MAX, vb + gv/GRAD_MAX);
            dataPack.idxVec.push_back(vs*img1.cols + us);
        }
    }
    // TODO check the reconstruction and discard bad points
    vector<Vector3d> cloudGradPointVec;
    depthMap.reconstruct(imageGradPointVec, cloudGradPointVec);
    depthMap.reconstruct(imagePointVec, dataPack.cloud);
    for (int i = 0; i < dataPack.cloud.size(); i++)
    {
        dataPack.gradientVec.push_back( (cloudGradPointVec[i] - dataPack.cloud[i]).normalized() );
        
    }
    return dataPack;
}

void ScalePhotometric::computePose(const Mat32f & img2, Transformation<double> & T12)
{
    if (verbosity > 0) 
    {
        cout << "ScalePhotometric::computePose" << endl;
    }
    scaleSpace2.generate(img2);
    for (int scaleIdx = scaleSpace1.size() - 1; scaleIdx >= 2; scaleIdx--)
    {
        computePose(scaleIdx, T12);
//        break;
    }
}

void ScalePhotometric::computePose(int scaleIdx, Transformation<double> & T12)
{
    if (verbosity > 1) 
    {
        cout << "ScalePhotometric::computePose with scaleIdx = " << scaleIdx << endl;
    }
    PhotometricPack dataPack = initPhotometricData(scaleIdx);
    scaleSpace2.setActiveScale(scaleIdx);
    const Mat32f & img2 = scaleSpace2.get();
    array<double, 6> pose = T12.toArray();
    Problem problem;
    PhotometricCostFunction * costFunction = new PhotometricCostFunction(camPtr2, dataPack,
                                                        img2, scaleSpace2.getActiveScale());
    problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), pose.data());
    
    
    //run the solver
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 150;
    if (verbosity > 2) options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    if (verbosity > 2) cout << summary.FullReport() << endl;
    else if (verbosity > 1) cout << summary.BriefReport() << endl;
    T12 = Transformation<double>(pose.data());
}

// Optimization using the auto differentiation
void ScalePhotometric::computePoseAuto(int scaleIdx, Transformation<double> & T12)
{
    if (verbosity > 1) 
    {
        cout << "ScalePhotometric::computePoseAuto with scaleIdx = " << scaleIdx << endl;
    }
    PhotometricPack dataPack = initPhotometricData(scaleIdx);
    scaleSpace2.setActiveScale(scaleIdx);
//    const Mat32f & img1 = scaleSpace1.get();
    const Mat32f & img2 = scaleSpace2.get();
    array<double, 6> pose = T12.toArray();
    typedef DynamicAutoDiffCostFunction<PhotometricError<EnhancedProjector>> photometricCF;
    PhotometricError<EnhancedProjector> * errorFunctor;
    errorFunctor = new PhotometricError<EnhancedProjector>(camPtr2->params, dataPack,
                             img2, scaleSpace2.getActiveScale());
    photometricCF * costFunctionAuto = new photometricCF(errorFunctor);
    costFunctionAuto->AddParameterBlock(6);
    costFunctionAuto->SetNumResiduals(dataPack.cloud.size());
    
    
    
    Problem problemAuto;
    problemAuto.AddResidualBlock(costFunctionAuto, new SoftLOneLoss(1), pose.data());
 
    //run the solver
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 150;
    if (verbosity > 2) options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problemAuto, &summary);
    if (verbosity > 2) cout << summary.FullReport() << endl;
    else if (verbosity > 1) cout << summary.BriefReport() << endl;
    T12 = Transformation<double>(pose.data());
}

//bool PhotometricLocalization::wrapImage(const Mat32f & img2, Mat32f & imgOut,
//        Transformation<double> & T12)
//{
//    imgOut.create(cam1.height, cam1.width);
//    imgOut.setTo(0);
//    
//    vector<Vector3d> cloud2;
//    T12.inverseTransform(dataPack.cloud, cloud2);
//    
//    vector<Vector2d> pointVec2;
//    cam2.projectPointCloud(cloud2, pointVec2);
//    
//    for (int i = 0; i < pointVec2.size(); i++)
//    {
//        int u = round(pointVec2[i][0]), v = round(pointVec2[i][1]);
//        if (u < 0 or u >= img2.cols or v < 0 or v > img2.rows) continue;
//        imgOut(dataPack.idxVec[i]) = img2(v, u);
//    }
//    return true;
//}

//const int step = 27;

//bool PhotometricLocalization::initCloud(const Mat32f & img1, const Mat32f & dist)
//{
//    vector<Vector2d> imagePointVec;
//    vector<double> distVec;
//    for (int v = 0; v < img1.rows; v += step)
//    {
//        for (int u = 0; u < img1.cols; u += step)
//        {
//            //TODO change it
//            int ud = params.uSmall(u);
//            int vd = params.vSmall(v);
//            if (ud < 0 or ud >= dist.cols or vd < 0 or vd >= dist.rows) continue;
//            if (dist(vd, ud) < 0.01) continue;
////            cout << u << " " << v <<" " << ud << " " << vd << " " << dist(vd, ud) << endl;
//            dataPack.colorVec.push_back(img1(v, u));
//            imagePointVec.emplace_back(u, v);
//            distVec.push_back(dist(vd, ud));
//            dataPack.idxVec.push_back(v*img1.cols + u);
//        }
//    }
////    cout << indexVec.size() << endl;
//    // TODO check the reconstruction and discard bad points
//    cam1.reconstructPointCloud(imagePointVec, dataPack.cloud);
//    for (int i = 0; i < dataPack.cloud.size(); i++)
//    {
//        dataPack.cloud[i] = dataPack.cloud[i] * (distVec[i] / dataPack.cloud[i].norm());
////        cout << cloud[i].transpose() << endl;
//    }
//    return true;
//}

