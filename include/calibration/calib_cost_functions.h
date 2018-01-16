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
    
#pragma once

#include "eigen.h"
#include "ceres.h"

#include "projection/generic_camera.h"

enum TransformationStatus {TRANSFORM_DIRECT, TRANSFORM_INVERSE};

struct GenericProjectionJac : ceres::CostFunction
{
    GenericProjectionJac(const Vector2dVec & proj, const Vector3dVec & grid,
            const ICamera * const camera,
            const vector<TransformationStatus> & transformStatusVec) : 
    _proj(proj),
    _grid(grid),
    _camera(camera->clone()),
    _transformStatusVec(transformStatusVec)
    {
//        int count = 0;
        mutable_parameter_block_sizes()->push_back(_camera->numParams()); //FIXME
        for (auto & x : _transformStatusVec)
        {
            mutable_parameter_block_sizes()->push_back(6);
        }
//        assert(_constTransformVec.size() == count);
//        mutable_parameter_block_sizes()->push_back(_camera->numParams()); //FIXME
        set_num_residuals(_proj.size()*2);
    }

    virtual ~GenericProjectionJac() 
    {
        delete _camera;
    }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    bool _variableIntrinsics;
    ICamera * _camera;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
    const vector<TransformationStatus> _transformStatusVec;
//    const vector<Transformation<double>> _constTransformVec;
};

struct OdometryPrior : ceres::SizedCostFunction<6, 6, 6>
{
    OdometryPrior(const double errV, const double errW, const double lambda,
        const Transformation<double> xi1,
        const Transformation<double> xi2);

    virtual ~OdometryPrior() { }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    Transf _zetaPrior;
    Matrix6d _A;
};

struct TransformationPrior : ceres::SizedCostFunction<6, 6>
{
    TransformationPrior(const double * const stiffness, const double * const xi):
        _xiPrior(xi),
        _A(Matrix6d::Zero()),
        _R(_xiPrior.rotMat())
    {
        for (int i = 0; i < 6; i++)
        {
            _A(i, i) = stiffness[i];
        }
        Matrix3d M = interOmegaRot(_xiPrior.rot());
        _A.bottomRightCorner<3, 3>() = _A.bottomRightCorner<3, 3>() * M;
        
    }

    virtual ~TransformationPrior() { }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const;
    
    Transf _xiPrior;
    Matrix6drm _A;
    Matrix3d _R;
};

