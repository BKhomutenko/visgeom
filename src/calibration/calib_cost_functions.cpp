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
    
#include "calibration/calib_cost_functions.h"

#include "eigen.h"
#include "io.h"

#include "geometry/geometry.h"
#include "projection/jacobian.h"
#include "projection/generic_camera.h"


bool GenericProjectionJac::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //compute the transformation chain
    Transformation<double> xiAcc;
    for (int paramIdx = 0; paramIdx < _transformStatusVec.size(); paramIdx++)
    {
        auto const status = _transformStatusVec[paramIdx];
        if (status == TRANSFORM_DIRECT)
        {
            Transformation<double> xi(params[paramIdx]);
            xiAcc = xiAcc.compose(xi);
        }
        else if (status == TRANSFORM_INVERSE)
        {
            Transformation<double> xi(params[paramIdx]);
            xiAcc = xiAcc.composeInverse(xi);
        }
    }
    
    //compute the points in the camera frame
    Vector3dVec pointCamVec;
    xiAcc.transform(_grid, pointCamVec);
    
    //get the intrinsic parameters
    const int cameraIntrinsicIdx = _transformStatusVec.size();
    _camera->setParameters(params[cameraIntrinsicIdx]);
    
    //compute the reprojection error
    for (int i = 0; i < pointCamVec.size(); i++)
    {
        Vector2d modProj;
        if (_camera->projectPoint(pointCamVec[i], modProj)) 
        {
            Vector2d diff = modProj - _proj[i];
            residual[2*i] = diff[0];
            residual[2*i + 1] = diff[1];
        }
        else
        {
            residual[2*i] = DOUBLE_BIG;
            residual[2*i + 1] = DOUBLE_BIG;
        }
    }
    
    if (jacobian != NULL)
    {
        //compute the transformation chain
        Transformation<double> xiAcc;
        for (int paramIdx = 0; paramIdx < _transformStatusVec.size(); paramIdx++)
        {
            auto const status = _transformStatusVec[paramIdx];
            Transformation<double> xi23(params[paramIdx]);
            Transformation<double> xi13;
            if (status == TRANSFORM_DIRECT)
            {
                xiAcc = xiAcc.compose(xi23);
                xi13 = xiAcc;
            }
            else if (status == TRANSFORM_INVERSE)
            {
                xi13 = xiAcc;
                xiAcc = xiAcc.composeInverse(xi23);
            }
            if (jacobian[paramIdx] != NULL)
            {
                InterJacobian jacobianCalculator(_camera, xi13, xi23, status == TRANSFORM_INVERSE);
                
                for (int i = 0; i < pointCamVec.size(); i++)
                {
                    jacobianCalculator.dpdxi(pointCamVec[i], jacobian[paramIdx] + i*12,
                                                            jacobian[paramIdx] + i*12 + 6);
                }
            }
        }
        
        if (jacobian[cameraIntrinsicIdx] != NULL)
        {
            
            for (int i = 0; i < pointCamVec.size(); i++)
            {
                _camera->intrinsicJacobian(pointCamVec[i],
                            jacobian[cameraIntrinsicIdx] + i * 2 * _camera->numParams(),
                            jacobian[cameraIntrinsicIdx] + (i * 2 + 1) * _camera->numParams());
            }
        }
    }
    return true;
}

OdometryPrior::OdometryPrior(const double errV, const double errW, const double lambda,
        const Transformation<double> xi1,
        const Transformation<double> xi2) : 
    _zetaPrior(xi1.inverseCompose(xi2)),
    _A(Matrix6d::Zero())
{
    
    const double MIN_SIGMA_V = 0.01;
    const double MIN_SIGMA_W = 0.01;
    
    const double MIN_DELTA = 0.01;
    const double MIN_L = 0.01;
    
    const double delta = max(_zetaPrior.rot().norm(), MIN_DELTA);
    const double l = max(_zetaPrior.trans().norm(), MIN_L);
    
    const double delta2 = delta / 2.;
    const double l2 = l / 2.;
    
    const double s = sin(delta2);
    const double c = cos(delta2);
    
    
    Matrixd<3, 2> dfdu;
    dfdu <<     c,      l2 * s,  
               -s,      l2 * c, 
                0,           1;
    
    Matrix2d Cu;
    Cu <<   errV * errV * l * l,              0,
                      0,    errW * errW * delta * delta; 
    Cu(0, 0) = max(Cu(0, 0), MIN_SIGMA_V*MIN_SIGMA_V);
    Cu(1, 1) = max(Cu(1, 1), MIN_SIGMA_W*MIN_SIGMA_W);
    
    Matrix3d Cx = dfdu * Cu * dfdu.transpose() + (lambda * lambda) * Matrix3d::Identity();
    Matrix3d CxInv = Cx.inverse();
    Eigen::LLT<Matrix3d> lltOfCxInv(CxInv);
    Matrix3d U = lltOfCxInv.matrixU();
    _A.topLeftCorner<2, 2>() = U.topLeftCorner<2, 2>();
    _A.topRightCorner<2, 1>() = U.topRightCorner<2, 1>();
    _A(2, 2) = 1. / lambda;
    
    Matrix3d B = Matrix3d::Zero();
    B(0, 0) = 1. / lambda;
    B(1, 1) = 1. / lambda;
    B(2, 2) = U(2, 2);
    
    _A.bottomRightCorner<3, 3>() = B /* interOmegaRot(_zetaPrior.rot())*/; //FIXME
    
//    _A = Matrix6d::Identity();
}

//TODO check whether the interaction matrix calculation is redundand here
//TODO replace (dxi - prior) by dxi.inverseCompose(prior)
bool OdometryPrior::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    Transformation<double> xi1(params[0]);
    Transformation<double> xi2(params[1]);
    Transformation<double> zeta = xi1.inverseCompose(xi2);
    
    
    Map<Vector6d> res(residual);
    Vector6d err;
    _zetaPrior.inverseCompose(zeta).toArray(err.data());
    res = _A * err;
    
    if (jacobian != NULL)
    {
        //first transformation
        if (jacobian[0] != NULL)
        {
            Matrix3d R10 = xi1.rotMatInv();
            Matrix3d R10_Mwr1 = R10 * interOmegaRot(xi1.rot());
            Matrix6d J1;
            J1 << R10, Matrix3d::Zero(), Matrix3d::Zero(), R10_Mwr1;
            Matrix6d TT = _zetaPrior.screwTransfInv();
            Map<Matrix6drm> jac(jacobian[0]);
            jac = -_A * TT * J1;
        }
        
        //second transformation
        if (jacobian[1] != NULL)
        {
            Matrix3d R20 = xi2.rotMatInv();
            Matrix3d R20_Mwr2 = R20 * interOmegaRot(xi2.rot());
            Matrix6d J2;
            J2 << R20, Matrix3d::Zero(), Matrix3d::Zero(), R20_Mwr2;
            Map<Matrix6drm> jac(jacobian[1]);
            jac = _A * J2;
        }
    }
    return true;    
}

bool TransformationPrior::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    Map<Vector6d> res(residual);
    Vector6d err;
    _xiPrior.inverseCompose( Transf(params[0]) ).toArray(err.data()) ;
    err.head<3>() = _R * err.head<3>();
    err.tail<3>() = _R * err.tail<3>();
    res = _A * err;
    if (jacobian != NULL and jacobian[0] != NULL)
    {
        copy(_A.data(), _A.data() + 36, jacobian[0]);
    }
    return true;
}
    
