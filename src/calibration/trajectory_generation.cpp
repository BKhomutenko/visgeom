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

#include "geometry/geometry.h"
    
TrajectoryQuality::TrajectoryQuality(ITrajectory * traj, Transf xiCam,
        const Matrix6d & CovVis, const Matrix6d & CovPrior) : 
    _traj(traj),
    _xiCam(xiCam),
    _infoVis(CovVis.inverse())
{
    Matrix6d infoPrior = CovPrior.inverse();
    Matrix3d M = interOmegaRot(xiCam.rot());
    Matrix3d Mt = M.transpose();
    
    _hessPrior.topLeftCorner<3, 3>() = infoPrior.topLeftCorner<3, 3>();  
    _hessPrior.topRightCorner<3, 3>() = infoPrior.topRightCorner<3, 3>() * M;
    _hessPrior.bottomLeftCorner<3, 3>() = Mt * infoPrior.bottomLeftCorner<3, 3>();
    _hessPrior.bottomRightCorner<3, 3>() = Mt * infoPrior.bottomRightCorner<3, 3>() * M;
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
    vector<Transf> xiOdomVec = _traj->compute(params);
    Matrix6d H = _hessPrior;
    for (auto & xi : xiOdomVec)
    {
        Transf xiOrigCam = xi.compose(_xiCam);
        Transf xiVis = _xiCam.inverseCompose(xiOrigCam);
        Matrix3d RcamOrig = xiOrigCam.rotMatInv();
        Matrix3d RorigCam = RcamOrig.transpose();
        Matrix6d infoOrig;
        infoOrig.topLeftCorner<3, 3>() = RorigCam * _infoVis.topLeftCorner<3, 3>() * RcamOrig;
        infoOrig.topRightCorner<3, 3>() = RorigCam * _infoVis.topRightCorner<3, 3>() * RcamOrig;
        infoOrig.bottomLeftCorner<3, 3>() = RorigCam * _infoVis.bottomLeftCorner<3, 3>() * RcamOrig;
        infoOrig.bottomRightCorner<3, 3>() = RorigCam * _infoVis.bottomRightCorner<3, 3>() * RcamOrig;
        Matrix6d jac = dxi1xi2dxi1(_xiCam, xiVis) - dxi1xi2dxi2(xi, _xiCam);
        H += jac.transpose() * infoOrig * jac;
    }
    JacobiSVD<Matrix6d> svd(H);
    double res = 0;
    for (int i = 0; i < 6; i++)
    {
        res -= log(svd.singularValues()(i));
    }
    return res;
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
