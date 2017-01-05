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
    _dxiPrior(),
    _A(Matrix6d::Zero())
{
    auto xi12 = xi1.inverseCompose(xi2);
    
    const double delta = xi12.rot().norm();
    const double l = xi12.trans().norm();
    
    const double delta2 = delta / 2.;
    const double l2 = l / 2.;
    
    const double s = sin(delta2);
    const double c = cos(delta2);
    
    xi12.toArray(_dxiPrior.data());
    
    Matrix32d dfdu;
    dfdu <<     c,     -l2 * s,  
                s,      l2 * c, 
                0,           1;
    
    Matrix2d Cu;
    Cu <<   errV * errV * l * l,              0,
                      0,    errW * errW * delta * delta; 
   
    Matrix3d Cx = dfdu * Cu * dfdu.transpose() + (lambda * lambda) * Matrix3d::Identity();
    Matrix3d CxInv = Cx.inverse();
    Eigen::LLT<Matrix3d> lltOfCxInv(CxInv); // compute the Cholesky decomposition of A
    Matrix3d U = lltOfCxInv.matrixU();
    _A.topLeftCorner<2, 2>() = U.topLeftCorner<2, 2>();
    _A.topRightCorner<2, 1>() = U.topRightCorner<2, 1>();
    _A(2, 2) = 1. / lambda;
    
    Matrix3d B = Matrix3d::Zero();
    B(0, 0) = 1. / lambda;
    B(1, 1) = 1. / lambda;
    B(2, 2) = U(2, 2);
    
    _A.bottomRightCorner<3, 3>() = B * interOmegaRot(xi12.rot());
}

//TODO check whether the interaction matrix calculation is redundand here
bool OdometryPrior::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    Transformation<double> xi1(params[0]);
    Transformation<double> xi2(params[1]);
    Transformation<double> xi12 = xi1.inverseCompose(xi2);
    
    Vector6d dxi(xi12.toArray().data());
    
    Map<Vector6d> res(residual);
    res = _A * (dxi - _dxiPrior);
    
    if (jacobian != NULL)
    {
        //common variables
        Matrix3d R10 = xi1.rotMatInv();
        Matrix3d Mrw = interRotOmega(xi12.rot());
        
        //first transformation
        if (jacobian[0] != NULL)
        {
            Matrix3d R10_Mwr1 = R10 * interOmegaRot(xi1.rot());
            Matrix3d t12skew = hat(xi12.trans());
            Matrix6d J1;
            J1 << -R10, t12skew * R10_Mwr1, Matrix3d::Zero(), -Mrw * R10_Mwr1;
            Map<Matrix6drm> jac(jacobian[0]);
            jac = _A * J1;
        }
        
        //second transformation
        if (jacobian[1] != NULL)
        {
            Matrix3d Mwr2 = interOmegaRot(xi2.rot());
            Matrix6d J2;
            J2 << R10, Matrix3d::Zero(), Matrix3d::Zero(), Mrw * R10 * Mwr2;
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
    res = _A * (Map<const Vector6d>(params[0]) - _xiPrior);
    if (jacobian != NULL and jacobian[0] != NULL)
    {
        copy(_A.data(), _A.data() + 36, jacobian[0]);
    }
    return true;
}
    
