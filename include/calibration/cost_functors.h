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

enum TransformationStatus {TRANFORM_DIRECT, TRANSFORM_INVERSE, TRANSFORM_CONSTANT};

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
    _constIntrinsicVec(intrinsicVec)  
    {
        assert(_constIntrinsicVec.size() == Projector<double>().countIntrinsic());
    }
    
    GenericProjection(const Vector2dVec & proj, const Vector3dVec & grid,
            const vector<TransformationStatus> & transformStatusVec,
            const vector<Transformation<double>> & constTransformVec) : 
    _proj(proj),
    _grid(grid),
    _transformStatusVec(transformStatusVec),
    _constTransformVec(constTransformVec)  {}
    
    GenericProjection(const Vector2dVec & proj, const Vector3dVec & grid,
            const vector<TransformationStatus> & transformStatusVec) : 
    _proj(proj),
    _grid(grid),
    _transformStatusVec(transformStatusVec)  {}
         
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
            if (status == TRANFORM_DIRECT)
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
        Vector3Vec<T> transformedPoints;
        transformedPoints.reserve(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            transformedPoints.emplace_back(_grid[i].template cast<T>());
        }
        xiAcc.transform(transformedPoints, transformedPoints);
        
        //get the intrinsic parameters
        Projector<T> projector;
        vector<T> intrinsicVec;
        if (_constIntrinsicVec.size() == 0)
        {
            intrinsicVec.resize(projector.countIntrinsic());
            copy(params[paramIdx], params[paramIdx] + projector.countIntrinsic(), intrinsicVec.begin());
        }
        else
        {
            for (auto & x : _constIntrinsicVec)
            {
                intrinsicVec.emplace_back(x);
            }
        }
        
        //compute the reprojection error
        for (int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(intrinsicVec.data(), transformedPoints[i].data(), modProj.data())) 
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
    const vector<double> _constIntrinsicVec;
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
        Vector3Vec<T> transformedPoints(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            transformedPoints[i] = _grid[i].template cast<T>();
        }
        TbaseGrid.transform(transformedPoints, transformedPoints);
        
        Projector<T> projector;
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(params[0], transformedPoints[i].data(), modProj.data())) 
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
        Vector3Vec<T> transformedPoints(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            transformedPoints[i] = _grid[i].template cast<T>();
        }
        TbaseGrid.transform(transformedPoints, transformedPoints);

        std::vector<T> camParamsT(_camParams.begin(), _camParams.end());
        Projector<T> projector;
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(camParamsT.data(), transformedPoints[i].data(), modProj.data()))
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

