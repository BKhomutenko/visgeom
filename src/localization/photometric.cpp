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
#include "ceres.h"

#include "projection/eucm.h"
#include "reconstruction/mh_pack.h"

void ScalePhotometric::setBaseImage(const Mat8u & img1)
{
    if (verbosity > 0) cout << "ScalePhotometric::computeBaseScaleSpace" << endl;
    scaleSpace1.generate(img1);
}

void ScalePhotometric::setTargetImage(const Mat8u & img2)
{
    if (verbosity > 0) cout << "ScalePhotometric::computeTargetScaleSpace" << endl;
    scaleSpace2.generate(img2);
}

void ScalePhotometric::setMotionPriorStatus(const bool val)
{
    useMotionPrior = val;
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
    vector<double> valVec;
    vector<int> packIdxVec;
    vector<Vector2d> imagePointVec;
    if (verbosity > 3) cout << "    scaled image size : " << img1.size() << endl;
    for (int vs = 0; vs < img1.rows; vs++)
    {
        for (int us = 0; us < img1.cols; us++)
        {
            double gu = gradU1(vs, us);
            double gv = gradV1(vs, us);
            if (gu*gu + gv*gv < GRAD_THRESH) continue; 
            
            //TODO change the threshold depending on the image size or use adaptive random sampling
            // to speed up the computation
            int ub = scaleSpace1.uConv(us);
            int vb = scaleSpace1.vConv(vs);
            if (depthMap.nearest(ub, vb) > DIST_MAX) continue;
            if (verbosity > 4) cout << "    " << vs << " " << us << endl;
            valVec.push_back(img1(vs, us));
            imagePointVec.emplace_back(ub, vb);
            packIdxVec.push_back(vs*img1.cols + us);
        }
    }
    vector<int> reconstIdxVec;
    depthMap.reconstruct(imagePointVec, reconstIdxVec, dataPack.cloud);
    _xiBaseCam.transform(dataPack.cloud, dataPack.cloud);
    for (auto & idx : reconstIdxVec)
    {
        dataPack.idxVec.push_back(packIdxVec[idx]);
        dataPack.valVec.push_back(valVec[idx]);
    }
    if (verbosity > 3) cout << "    datapack size : " << reconstIdxVec.size() << endl;
    return dataPack;
}

Transf ScalePhotometric::computePose(const Transf & T12)
{
    if (verbosity > 0) 
    {
        cout << "ScalePhotometric::computePose" << endl;
    }
    //TODO set the optimization depth with the parameters   v
    Transf xi = T12;
    _xiPrior = T12;
    for (int scaleIdx = scaleSpace1.size() - 1; scaleIdx >= 0; scaleIdx--)
    {
        computePose(scaleIdx, xi);
    }
    return xi;
}

void ScalePhotometric::computePose(int scaleIdx, Transf & T12)
{
    if (verbosity > 1) 
    {
        cout << "ScalePhotometric::computePose with scaleIdx = " << scaleIdx << endl;
    }
    PhotometricPack dataPack = initPhotometricData(scaleIdx);
    scaleSpace2.setActiveScale(scaleIdx);
    array<double, 6> pose = T12.toArray();
    Problem problem;
    PhotometricCostFunction * costFunction = new PhotometricCostFunction(camPtr2, _xiBaseCam, dataPack,
                                            scaleSpace2.get(), scaleSpace2.getActiveScale());
    
    const double LOSS_THRESH = 10;
//    RobustLoss myLoss(LOSS_THRESH);
//    CauchyLoss myLoss(LOSS_THRESH);
//    for (int i = 0; i < 100; i+=5)
//    {
//        double out[3];
//        myLoss.Evaluate(i*i, out);
//        cout << setw(15) << i      << setw(15) << out[0]
//             << setw(15) << out[1] << setw(15) << out[2] << endl;
//    }
    
//    problem.AddResidualBlock(costFunction, new RobustLoss(LOSS_THRESH), pose.data());
//    problem.AddResidualBlock(costFunction, new CauchyLoss(LOSS_THRESH), pose.data());
    problem.AddResidualBlock(costFunction, NULL, pose.data());    
    
    //add an odometry prior
    if (useMotionPrior and 0) //FIXME experimental
    {
        //A proper nose model on the depth map localization must be applied
        OdometryPrior * odometryCost = new OdometryPrior(0.03, 0.03, 0.03, 0.03, _xiPrior);
        problem.AddResidualBlock(odometryCost, NULL, pose.data());
    }
    
    //run the solver
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 150;
    if (verbosity > 2) options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    if (verbosity > 2) cout << summary.FullReport() << endl;
    else if (verbosity > 1) cout << summary.BriefReport() << endl;
    T12 = Transf(pose.data());
}

Transf ScalePhotometric::computePoseMI(const Transf & T12)
{
    if (verbosity > 0) 
    {
        cout << "ScalePhotometric::computePoseMI" << endl;
    }
    Transf xi = T12;
    //TODO set the optimization depth with the parameters   v
    for (int scaleIdx = scaleSpace1.size() - 1; scaleIdx >= 0; scaleIdx--)
    {
        computePoseMI(scaleIdx, xi);
    }
    return xi;
}

//TODO put to a separate file   
void saveSurface(string fileName, FirstOrderFunction * func, 
        int idx1, int idx2, double step, int Nsteps, double* params)
{
    vector<double> paramVec(params, params + func->NumParameters());
    ofstream surfaceFile;
    surfaceFile.open(fileName);
    for (int i = -Nsteps; i <= Nsteps; i++)
    {
        for (int j = -Nsteps; j <= Nsteps; j++)
        {
            double f;
            paramVec[idx1] = params[idx1] + i*step;
            paramVec[idx2] = params[idx2] + j*step;
            func->Evaluate(paramVec.data(), &f, NULL);
            surfaceFile << setw(15) << f;
        }
        surfaceFile << endl;
    }
    surfaceFile.close();
}

array<double, 6> ScalePhotometric::covarianceEigenValues(const int scaleIdx,
        const Transf T12, bool baseValues)
{
    PhotometricPack dataPack = initPhotometricData(scaleIdx);
    PhotometricCostFunction * costFunction;
    if (baseValues)
    {
        //FIXME must be camPtr1
        costFunction = new PhotometricCostFunction(camPtr2, _xiBaseCam, dataPack,
                                scaleSpace1.get(), scaleSpace1.getActiveScale());
    }
    else
    {
        scaleSpace2.setActiveScale(scaleIdx);
        costFunction = new PhotometricCostFunction(camPtr2, _xiBaseCam, dataPack,
                                scaleSpace2.get(), scaleSpace2.getActiveScale());
    }
    cout << "Cost function is created" << endl;
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> J(dataPack.valVec.size(), 6);
    
    //TODO fix this mess
    vector<double> residual(dataPack.valVec.size());
    array<double, 6> pose = T12.toArray();
    double * jacPtr = J.data();
    double * paramPtr = pose.data();
    costFunction->Evaluate(&paramPtr, residual.data(), &jacPtr);
    
    vector<int> errHist(255);
    for(double err : residual)
    {
        int idx = min( 255, max( 0, int(round(abs(err))) ) );
        errHist[idx]++;
    }
    
//    for (int h : errHist)
//    {
//        cout << h << endl;
//    }
    cout << "Jacobian is evaluated" << endl;
    MatrixXd JtJ = J.transpose() * J;
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(JtJ);
    const double * evPtr = es.eigenvalues().data();
    array<double, 6> res;
    copy(evPtr, evPtr + 6, res.data());
    delete costFunction;
    return res;    
}

void ScalePhotometric::computePoseMI(int scaleIdx, Transf & T12)
{
    if (verbosity > 1) 
    {
        cout << "ScalePhotometric::computePoseMI with scaleIdx = " << scaleIdx << endl;
    }
    PhotometricPack dataPack = initPhotometricData(scaleIdx);
    scaleSpace2.setActiveScale(scaleIdx);
    array<double, 6> pose = T12.toArray();
    MutualInformation * costFunction = new MutualInformation(camPtr2, dataPack, _xiBaseCam,
                                scaleSpace2.get(), scaleSpace2.getActiveScale(), 8, 255);
    
//    if (useMotionPrior) //FIXME experimental
//    {
//        //A proper nose model on the depth map localization must be applied
//        OdometryPrior * odometryCost = new OdometryPrior(0.03, 0.03, 0.01, 0.03, _xiPrior);
//        problem.AddResidualBlock(odometryCost, NULL, pose.data());
//    }
    
    //TODO put to a separate file                                                    
//    array<double, 6> grad;
//    double val;
//    array<double, 6> val1;
//    array<double, 6> val2;
//    double eps = 1e-5;
//    costFunction->Evaluate(pose.data(), &val, grad.data());
//    for (int i = 0; i < 6; i++)
//    {
//        pose[i] += eps;
//        costFunction->Evaluate(pose.data(), &val2[i], NULL);
//        pose[i] -= 2*eps;
//        costFunction->Evaluate(pose.data(), &val1[i], NULL);
//        pose[i] += eps;
//    }
//    
//    for (int i = 0; i < 6; i++)
//    {
//        cout << grad[i] << "    " << (val2[i] - val1[i]) / 2 / eps << endl;
//    }
    
//    return;
    GradientProblem problem(costFunction);
    
    if (verbosity > 2) cout << "    Problem created" << endl;
    //run the solver
    GradientProblemSolver::Options options;
    options.line_search_direction_type = ceres::BFGS;
//    options.use_approximate_eigenvalue_bfgs_scaling = true;
//    options.line_search_interpolation_type = ceres::CUBIC;
    options.function_tolerance = 1e-2;
    options.gradient_tolerance = 1e-3;
//    options.linear_solver_type = ceres::DENSE_QR;
//    options.max_num_iterations = 15;
    if (verbosity > 2) options.minimizer_progress_to_stdout = true;
    GradientProblemSolver::Summary summary;
    Solve(options, problem, pose.data(), &summary);
    if (verbosity > 2) cout << summary.FullReport() << endl;
    else if (verbosity > 1) cout << summary.BriefReport() << endl;
    T12 = Transf(pose.data());
//    cout << T12 << endl;
//    saveSurface("surf01.txt", costFunction, 2, 3, 0.0005, 50, pose.data());
}

void ScalePhotometric::wrapImage(const Mat8u & src, Mat8u & dst, const Transf T12) const
{   
    MHPack pack;
    dst.create(src.size());
    dst.setTo(0);
    for (int i = 0; i < src.rows; i++)
    {
        for (int j = 0; j < src.cols; j++)
        {
            pack.imagePointVec.emplace_back(j, i);
        }
    }
    //TODO rewrite the reconstruct
    depthMap.reconstruct(pack, QUERY_POINTS | INDEX_MAPPING);
    Transf xiCam12 = _xiBaseCam.inverseCompose(T12).compose(_xiBaseCam);
    xiCam12.inverseTransform(pack.cloud, pack.cloud);
    Vector2dVec ptVec;
    depthMap.project(pack.cloud, ptVec);
    for (int i = 0; i < pack.cloud.size(); i++)
    {
        int u = round(ptVec[i][0]);
        int v = round(ptVec[i][1]);
        if (u < 0 or u >= src.cols) continue;
        if (v < 0 or v >= src.rows) continue;
        int v1 = pack.idxMapVec[i] / src.cols;
        int u1 = pack.idxMapVec[i] % src.cols; // FIXME must be the other way around
        dst(v1, u1) = src(v, u);
    }
}
