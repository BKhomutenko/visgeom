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

#include "utils/jacobian.h"
#include "camera/generic_camera.h"
#include "eigen.h"
#include "io.h"

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
        int count = 0;
        for (auto & x : _transformStatusVec)
        {
            mutable_parameter_block_sizes()->push_back(6);
        }
        assert(_constTransformVec.size() == count);
        mutable_parameter_block_sizes()->push_back(_camera->numParams());
        set_num_residuals(_proj.size()*2);
    }

    virtual ~GenericProjectionJac() 
    {
        delete _camera;
    }
    
    virtual bool Evaluate(double const * const * params,
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
                residual[2*i] = 0;
                residual[2*i + 1] = 0;
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
    
    bool _variableIntrinsics;
    ICamera * _camera;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
    const vector<TransformationStatus> _transformStatusVec;
    const vector<Transformation<double>> _constTransformVec;
};

template<template<typename> class Projector>
struct GridProjection 
{
    GridProjection(const Vector2dVec & proj, const Vector3dVec & grid)
    : _proj(proj), _grid(grid) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> TbaseGrid(params[1]);
        Vector3Vec<T> pointCamVec(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            pointCamVec[i] = _grid[i].template cast<T>();
        }
        TbaseGrid.transform(pointCamVec, pointCamVec);
        
        Projector<T> projector;
        for (unsigned int i = 0; i < pointCamVec.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(params[0], pointCamVec[i].data(), modProj.data())) 
            {
                Vector2<T> diff = _proj[i].template cast<T>() - modProj;
                residual[2*i] = diff[0];
                residual[2*i + 1] = diff[1];
            }
            else
            {
                residual[2*i] = T(0.);
                residual[2*i + 1] = T(0.);
            }
        }
        return true;
    }
    
    const Vector2dVec _proj;
    const Vector3dVec _grid;
};

    
template<template<typename> class Projector>
struct GridEstimate
{
    GridEstimate(const Vector2dVec & proj, const Vector3dVec & grid,
    const std::vector<double> & camParams) : _proj(proj), _grid(grid), _camParams(camParams) {}
            
    template <typename T>
    bool operator()(const T * const * params,
                    T* residual) const 
    {
        Transformation<T> TbaseGrid(params[0]);
        Vector3Vec<T> pointCamVec(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            pointCamVec[i] = _grid[i].template cast<T>();
        }
        TbaseGrid.transform(pointCamVec, pointCamVec);

        std::vector<T> camParamsT(_camParams.begin(), _camParams.end());
        Projector<T> projector;
        for (unsigned int i = 0; i < pointCamVec.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(camParamsT.data(), pointCamVec[i].data(), modProj.data()))
            {
                Vector2<T> diff = _proj[i].template cast<T>() - modProj;
                residual[2*i] = diff[0];
                residual[2*i + 1] = diff[1];
            }
            else
            {
                residual[2*i] = T(0.);
                residual[2*i + 1] = T(0.);
            }
        }
        return true;
    }
    
    const std::vector<double> _camParams;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
};

