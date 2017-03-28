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

#include "std.h"
#include "io.h"
#include "eigen.h"
#include "ceres.h"
#include "json.h"

#include "projection/generic_camera.h"
#include "geometry/geometry.h"

class ITrajectory
{
public:
    ITrajectory() {}
    virtual ~ITrajectory() {}
    virtual void compute(const double * params, vector<Transf> & trajVec, 
            vector<Matrix6d> & covVec) const = 0;
    virtual int paramSize() const = 0;
    
};

const double DIFF_EPS = 1e-8; //FIXME

struct TrajectoryQuality : FirstOrderFunction
{
    // the class takes the ownership of traj
    TrajectoryQuality(ITrajectory * traj, Transf xiCam,
            const Matrix6d & CovVis, const Matrix6d & CovPrior);
    
    virtual bool Evaluate(const double * params,
            double * residual, double * jacobian) const;
    
    double EvaluateCost(const double * params) const;
    
    virtual int NumParameters() const { return _traj->paramSize(); }
    virtual ~TrajectoryQuality()
    {
        delete _traj;
    }
    
    ITrajectory * _traj;
    Transf _xiCam;
    Matrix6d _covVis;  // invers of Covariance_vis
    Matrix6d _hessPrior;  // C_prior^-1   --  the regularization term
};

//Takes into account the field-of-view constraint and the visual localization quality
struct TrajectoryVisualQuality : FirstOrderFunction
{
    // the class takes the ownership of traj
    TrajectoryVisualQuality(const vector<ITrajectory*> & trajVec, const ICamera * camera,
            Transf xiCam, Transf xiBoard, const Vector3dVec & board,
            const Matrix6d & CovPrior, const Matrix2d & CovPt, 
            const int Nx, const int Ny, double kapaMax = 0.3);
    
    virtual bool Evaluate(const double * params,
            double * residual, double * jacobian) const;
    
    Matrix6d visualCov(const Transf & camPose) const;
    
    //TODO to think how to regularize for different cameras
    //This function is adapted for fisheye cameras
    double imageLimitsCost(const Transf & camPose) const;
    
    double curvatureCost(const Transf & xi1, const Transf & xi2) const;
    
    double EvaluateCost(const double * params) const;
    
    virtual int NumParameters() const { return _paramSize; }
    virtual ~TrajectoryVisualQuality()
    {
        for (auto & traj : _trajVec)
        {
            delete traj;   
        }
        delete _camera;
    }
    
    ICamera * _camera;
    vector<ITrajectory*> _trajVec;
    int _paramSize;
    Transf _xiCam, _xiBoard;
    Vector3dVec _board;
    int _Nx, _Ny;
    Matrix2d _ptStiffness;  // A^T * A = C_pt^-1
    Matrix6d _hessPrior;  // C_prior^-1   --  the regularization term
    double _kapaMax;
};

//TODO put lesewhere
// Transformation jacobian
Matrix6d dxi1xi2dxi1(const Transf & xi1, const Transf & xi2);

// Transformation jacobian
Matrix6d dxi1xi2dxi2(const Transf & xi1, const Transf & xi2);

