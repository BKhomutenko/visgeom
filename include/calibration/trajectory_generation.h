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
    Matrix6d _hessPrior;  // L^T * C_prior^-1 * L   --  the regularization term
};


//TODO put lesewhere
// Transformation jacobian
Matrix6d dxi1xi2dxi1(const Transf & xi1, const Transf & xi2);

// Transformation jacobian
Matrix6d dxi1xi2dxi2(const Transf & xi1, const Transf & xi2);

