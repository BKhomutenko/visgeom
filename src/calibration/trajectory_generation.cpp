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

#include "calibration/trajectory_generation.h"

#include "eigen.h"

#include "geometry/geometry.h"
  
/////////////////////////
/// TrajectoryQuality ///
/////////////////////////

TrajectoryQuality::TrajectoryQuality(ITrajectory * traj, Transf xiCam,
        const Matrix6d & CovVis, const Matrix6d & CovPrior) : 
    _traj(traj),
    _xiCam(xiCam),
    _covVis(CovVis),
    _hessPrior(CovPrior.inverse())
{
}
    
bool TrajectoryQuality::Evaluate(const double * params,
        double * residual, double * jacobian) const
{
    residual[0] = EvaluateCost(params);
    vector<double> paramsDiff(params, params + _traj->paramSize());
    for (int i = 0; i < paramsDiff.size(); i++)
    {
        //numeric differentiation
        double delta = max(abs(paramsDiff[i]) * DIFF_EPS, DIFF_EPS * DIFF_EPS);
        paramsDiff[i] += delta;
        double vPlus = EvaluateCost(paramsDiff.data());
        paramsDiff[i] -= 2*delta;
        double vMinus = EvaluateCost(paramsDiff.data());
        paramsDiff[i] += delta;
        jacobian[i] = (vPlus - vMinus) / (2 * delta);  
    }
    return true;
}


double TrajectoryQuality::EvaluateCost(const double * params) const
{
    vector<Transf> xiOdomVec;
    vector<Matrix6d> covOdomVec;
    _traj->compute(params, xiOdomVec, covOdomVec);
    Matrix6d H = _hessPrior;
    for (int i = 0; i < xiOdomVec.size(); i++)
    {
        Matrix6d LcamOdom = _xiCam.screwTransfInv();
        Matrix6d C = _covVis + LcamOdom * covOdomVec[i] * LcamOdom.transpose();
        Transf xiVis = _xiCam.inverseCompose(xiOdomVec[i]).compose(_xiCam);
        Matrix6d J = xiVis.screwTransfInv() - Matrix6d::Identity();
        // hessian is computed up to an orthonormal transformation
        // it does not change the rank but simplifies the calculus
        H += (J.transpose() * C.inverse() * J); 
    }
    JacobiSVD<Matrix6d> svd(H);
    double res = 0;
    for (int i = 0; i < 6; i++)
    {
        res -= log(svd.singularValues()(i));
    }
    return res;
}

///////////////////////////////
/// TrajectoryVisualQuality ///
///////////////////////////////

TrajectoryVisualQuality::TrajectoryVisualQuality(const vector<ITrajectory*> & trajVec, const ICamera * camera,
            Transf xiCam, Transf xiBoard, const Vector3dVec & board,
            const Matrix6d & CovPrior, const Matrix2d & CovPt, 
            const int Nx, const int Ny, double kapaMax) : 
    _trajVec(trajVec),
    _camera(camera->clone()),
    _xiCam(xiCam),
    _xiBoard(xiBoard),
    _hessPrior(CovPrior.inverse()),
    _board(board),
    _paramSize(0),
    _kapaMax(kapaMax),
    _Nx(Nx),
    _Ny(Ny)
{
    Eigen::LLT<Matrix2d> lltOfCovCornerInv(CovPt.inverse()); // compute the Cholesky decomposition of A
    _ptStiffness = lltOfCovCornerInv.matrixU();
    for (auto & t : _trajVec)
    {
        _paramSize += t->paramSize();
    }
}
    
bool TrajectoryVisualQuality::Evaluate(const double * params,
        double * residual, double * jacobian) const
{
    residual[0] = EvaluateCost(params);
    vector<double> paramsDiff(params, params + _paramSize);
    for (int i = 0; i < paramsDiff.size(); i++)
    {
        //numeric differentiation
        double delta = max(abs(paramsDiff[i]) * DIFF_EPS, DIFF_EPS * DIFF_EPS);
        double origValue = paramsDiff[i];
        paramsDiff[i] = origValue + delta;
        double vPlus = EvaluateCost(paramsDiff.data());
        paramsDiff[i] = origValue - delta;
        double vMinus = EvaluateCost(paramsDiff.data());
        paramsDiff[i] = origValue;
        jacobian[i] = (vPlus - vMinus) / (2 * delta);  
    }
    return true;
}

double TrajectoryVisualQuality::EvaluateCost(const double * params) const
{
    vector<Transf> xiOdomVec;
    vector<Matrix6d> covOdomVec;
    Matrix6d H = _hessPrior;
    double res = 0;
    for (int trajIdx = 0; trajIdx < _trajVec.size(); trajIdx++)
    {
        auto traj = _trajVec[trajIdx];
        traj->compute(params + trajIdx * traj->paramSize(), xiOdomVec, covOdomVec);
        
        for (int i = 1; i < xiOdomVec.size(); i++)
        {
            Matrix6d LcamOdom = _xiCam.screwTransfInv();
            Transf xiBaseCam = xiOdomVec[i].compose(_xiCam);
            Matrix6d C = visualCov(xiBaseCam) + LcamOdom * covOdomVec[i] * LcamOdom.transpose();
            Transf xiVis = _xiCam.inverseCompose(xiOdomVec[0].inverseCompose(xiBaseCam));
            Matrix6d J = xiVis.screwTransfInv() - Matrix6d::Identity();
            // hessian is computed up to an orthonormal transformation
            // it does not change the rank but simplifies the calculus
            H += (J.transpose() * C.inverse() * J); 
            res += imageLimitsCost(xiBaseCam);
            res += curvatureCost(xiOdomVec[i - 1], xiOdomVec[i]);
            //TODO add the distance to board as a constraint
    //        res += distanceCost();
        }
    }
    JacobiSVD<Matrix6d> svd(H);
    for (int i = 0; i < 6; i++)
    {
        res -= log(svd.singularValues()(i));
    }
    return res;
}

Matrix6d TrajectoryVisualQuality::visualCov(const Transf & camPose) const
{
    Transf xiCamBoard = camPose.inverseCompose(_xiBoard);
    Vector3dVec boardCamVec;
    xiCamBoard.transform(_board, boardCamVec);
    
    Matrix6d JtCJ = Matrix6d::Zero(), JtJ = Matrix6d::Zero();
    for (int i = 0; i < boardCamVec.size(); i++)
    {
        Matrix23drm dpdx;
        _camera->projectionJacobian(boardCamVec[i], dpdx.data(), dpdx.data() + 3);
        
        Matrixd<3, 6> dxdxi;
        dxdxi.topLeftCorner<3, 3>() = -Matrix3d::Identity();
        dxdxi.topRightCorner<3, 3>() = hat(boardCamVec[i]);
        Matrixd<2, 6> dpdxi = dpdx * dxdxi;
        JtJ += dpdxi.transpose() * dpdxi;
        JtCJ += dpdxi.transpose() * _ptStiffness * dpdxi;
    }
    Matrix6d JtCJinv = JtCJ.inverse();
    return JtCJinv.transpose() * JtJ * JtCJinv;
}

//image borders and the board size
double TrajectoryVisualQuality::imageLimitsCost(const Transf & camPose) const
{
    const double MARGIN = 100;
    const double MIN_DIAGONAL = 50;
    //FIXME temporary for particular data
    Vector2d center(_camera->getCenterU(), _camera->getCenterV());
//    Vector2d center(280, 215);
//    double r = min(center[1], _camera->height - center[1])  - MARGIN; // margin
    double r = (_camera->getCenterU() + _camera->getCenterV()) / 2;
    double res = 0;
    
    Transf xiCamBoard = camPose.inverseCompose(_xiBoard);
    Vector3dVec boardCamVec;
    xiCamBoard.transform(_board, boardCamVec);
    
    Vector2dVec boardProjectionVec;
    _camera->projectPointCloud(boardCamVec, boardProjectionVec);
    for (auto & pt : boardProjectionVec)
    {
        double distX = max(MARGIN - min(pt[0], _camera->width - pt[0]), 0.);
        double distY = max(MARGIN - min(pt[1], _camera->height - pt[1]), 0.);
        double distR = max((pt - center).norm() - r, 0.);
        res += pow(max(distX, max(distY, distR))/ 10., 2);
    }
    auto pt1 = boardProjectionVec[0];
    auto pt2 = boardProjectionVec.back();
    auto pt3 = boardProjectionVec[_Nx - 1];
    auto pt4 = boardProjectionVec[_Nx * (_Ny- 1)];
    double diag1 = (pt1 - pt2).norm();
    double diag2 = (pt3 - pt4).norm();
    double diag = min(diag1, diag2);
    res += pow(max(MIN_DIAGONAL - (pt1 - pt2).norm(), 0.), 2);
    res += pow(max(abs(diag1 - diag2) - MIN_DIAGONAL / 2, 0.), 2);
    return res;
}

double TrajectoryVisualQuality::curvatureCost(const Transf & xi1, const Transf & xi2) const
{
    const Transf zeta = xi1.inverseCompose(xi2);
    const double v = zeta.trans().norm();
    const double w = zeta.rot().norm();
    const double kapa = w / v;
    return pow(30 * max(0., kapa - _kapaMax), 2);
}

// Transformation kinematic jacobian
Matrix6d dxi1xi2dxi1(const Transf & xi1, const Transf & xi2)
{
    Matrix6d jac;
    jac.topLeftCorner<3, 3>() = Matrix3d::Identity();
    jac.bottomLeftCorner<3, 3>() = Matrix3d::Zero();
    Matrix3d Momega = interOmegaRot(xi1.rot());
    Transf xi1xi2 = xi1.compose(xi2);
    Vector3d t = xi1.trans() - xi1xi2.trans();
    jac.topRightCorner<3, 3>() = hat(t) * Momega;
    jac.bottomRightCorner<3, 3>() = Momega;
    return jac;
}

// Transformation kinematic jacobian
Matrix6d dxi1xi2dxi2(const Transf & xi1, const Transf & xi2)
{
    Matrix6d jac;
    Matrix3d R01 = xi1.rotMat();
    Transf xi1xi2 = xi1.compose(xi2);
    jac.topRightCorner<3, 3>() = Matrix3d::Zero();
    jac.topLeftCorner<3, 3>() = R01;
    jac.bottomLeftCorner<3, 3>() = Matrix3d::Zero();
    jac.bottomRightCorner<3, 3>() = R01 * interOmegaRot(xi2.rot());
    return jac;
}
