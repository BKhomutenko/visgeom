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

enum TransformationStatus {TRANSFORM_DIRECT, TRANSFORM_INVERSE, TRANSFORM_CONSTANT};

template<template<typename> class Projector>
struct GenericProjection
{
    GenericProjection(const Vector2dVec & proj, const Vector3dVec & grid,
            const vector<TransformationStatus> & transformStatusVec,
            const vector<Transformation<double>> & constTransformVec,
            const vector<double> & intrinsicVec) : 
    _proj(proj),
    _grid(grid),
    _transformStatusVec(transformStatusVec),
    _constTransformVec(constTransformVec),
    _constIntrinsicVec(intrinsicVec),
    _intrinsicCount(Projector<double>::INTRINSIC_COUNT)
    {
        assert(_constIntrinsicVec.size() == _intrinsicCount or 
               _constIntrinsicVec.size() == 0);
    }
    
    GenericProjection(const Vector2dVec & proj, const Vector3dVec & grid,
            const vector<TransformationStatus> & transformStatusVec,
            const vector<Transformation<double>> & constTransformVec) : 
    _proj(proj),
    _grid(grid),
    _transformStatusVec(transformStatusVec),
    _constTransformVec(constTransformVec),
    _intrinsicCount(Projector<double>::INTRINSIC_COUNT)  {}
    
    GenericProjection(const Vector2dVec & proj, const Vector3dVec & grid,
            const vector<TransformationStatus> & transformStatusVec) : 
    _proj(proj),
    _grid(grid),
    _transformStatusVec(transformStatusVec),
    _intrinsicCount(Projector<double>::INTRINSIC_COUNT)  {}
         
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        int paramIdx = 0;
        int constIdx = 0;
        
        //compute the transformation chain
        Transformation<T> xiAcc;
        for (TransformationStatus status : _transformStatusVec)
        {
            if (status == TRANSFORM_DIRECT)
            {
                Transformation<T> xi(params[paramIdx]);
                xiAcc = xiAcc.compose(xi);
                paramIdx++;
            }
            else if (status == TRANSFORM_INVERSE)
            {
                Transformation<T> xi(params[paramIdx]);
                xiAcc = xiAcc.composeInverse(xi);
                paramIdx++;
            }
            else if (status == TRANSFORM_CONSTANT)
            {
                Transformation<T> xi = _constTransformVec[constIdx].template cast<T>();
                xiAcc = xiAcc.compose(xi);
                constIdx++;
            }
        }
        
        //compute the points in the camera frame
        Vector3Vec<T> pointCamVec;
        pointCamVec.reserve(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            pointCamVec.emplace_back(_grid[i].template cast<T>());
        }
        xiAcc.transform(pointCamVec, pointCamVec);
        
        //get the intrinsic parameters
        Projector<T> projector;
        vector<T> intrinsicVec;
        if (_constIntrinsicVec.size() == 0)
        {
            intrinsicVec.resize(_intrinsicCount);
            copy(params[paramIdx], params[paramIdx] + _intrinsicCount, intrinsicVec.begin());
        }
        else
        {
            intrinsicVec.reserve(_intrinsicCount);
            for (auto & x : _constIntrinsicVec)
            {
                intrinsicVec.emplace_back(x);
            }
        }
        
        //compute the reprojection error
        for (int i = 0; i < pointCamVec.size(); i++)
        {
            Vector2<T> modProj;
            if (pointCamVec[i] != pointCamVec[i])
            {
                cout << i << endl;
            }
            if (projector(intrinsicVec.data(), pointCamVec[i].data(), modProj.data())) 
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
    
    const int _intrinsicCount;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
    const vector<double> _constIntrinsicVec;
    const vector<TransformationStatus> _transformStatusVec;
    const vector<Transformation<double>> _constTransformVec;
};

struct GenericProjectionJac : ceres::CostFunction
{
    GenericProjectionJac(const Vector2dVec & proj, const Vector3dVec & grid,
            const ICamera * const camera,
            const vector<TransformationStatus> & transformStatusVec,
            const vector<Transformation<double>> & constTransformVec = vector<Transformation<double>>(),
            bool variableIntrinsics = false) : 
    _proj(proj),
    _grid(grid),
    _camera(camera->clone()),
    _variableIntrinsics(variableIntrinsics),
    _transformStatusVec(transformStatusVec),
    _constTransformVec(constTransformVec)
    {
        int count = 0;
        for (auto & x : _transformStatusVec)
        {
            if (x == TRANSFORM_CONSTANT) count++;
            else mutable_parameter_block_sizes()->push_back(6);
        }
        assert(_constTransformVec.size() == count);
        if (variableIntrinsics) mutable_parameter_block_sizes()->push_back(_camera->numParams());
        set_num_residuals(_proj.size()*2);
    }

    virtual ~GenericProjectionJac() 
    {
        delete _camera;
    }
    
    virtual bool Evaluate(double const * const * params,
            double * residual, double ** jacobian) const
    {
        int paramIdx = 0;
        int constIdx = 0;
        
        //compute the transformation chain
        Transformation<double> xiAcc;
        for (TransformationStatus status : _transformStatusVec)
        {
            if (status == TRANSFORM_DIRECT)
            {
                Transformation<double> xi(params[paramIdx]);
                xiAcc = xiAcc.compose(xi);
                paramIdx++;
            }
            else if (status == TRANSFORM_INVERSE)
            {
                Transformation<double> xi(params[paramIdx]);
                xiAcc = xiAcc.composeInverse(xi);
                paramIdx++;
            }
            else if (status == TRANSFORM_CONSTANT)
            {
                xiAcc = xiAcc.compose(_constTransformVec[constIdx]);
                constIdx++;
            }
        }
        
        //compute the points in the camera frame
        Vector3dVec pointCamVec;
        xiAcc.transform(_grid, pointCamVec);
        
        //get the intrinsic parameters
        if (_variableIntrinsics)
        {
            _camera->setParameters(params[paramIdx]);
        }
        
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
                residual[2*i] = DOUBLE_INF;
                residual[2*i + 1] = DOUBLE_INF;
            }
        }
        
        if (jacobian != NULL)
        {
            int paramIdx = 0;
            int constIdx = 0;
            
            //compute the transformation chain
            Transformation<double> xiAcc;
            for (TransformationStatus status : _transformStatusVec)
            {
                if (status == TRANSFORM_CONSTANT)
                {
                    xiAcc = xiAcc.compose(_constTransformVec[constIdx]);
                    constIdx++;
                    continue;
                }
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
                paramIdx++;
            }
            
            if (_variableIntrinsics and jacobian[paramIdx] != NULL)
            {
                for (int i = 0; i < pointCamVec.size(); i++)
                {
                    _camera->intrinsicJacobian(pointCamVec[i],
                                    jacobian[paramIdx] + i * 2 * _camera->numParams(),
                                    jacobian[paramIdx] + (i * 2 + 1) * _camera->numParams());
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

