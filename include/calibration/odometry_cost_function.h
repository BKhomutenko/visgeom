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

/*
Cost functions for localization based on photometric data and mutual information
*/

#pragma once

#include "std.h"
#include "eigen.h"
#include "ocv.h"
#include "ceres.h"
#include "array"

#include "geometry/geometry.h"

//TODO refactor
struct OdometryCost : ceres::SizedCostFunction<6, 6>
{

public:
    
    OdometryCost(const double errV, const double errW, const double lambda,
                           const vector <Vector2d> & delta_q_vec, const double * const intrinsic_prior);

    virtual ~OdometryCost() { }
    
    virtual bool Evaluate(double const * const * params,
                          double * residual, double ** jacobian) const;

    void tf0n_jac_calc (const double * const odom_intrin,
                                  vector <Matrix3d> & jac_zeta_vec, 
                                  vector <Transf> & tf0_n_vec ) const;

    Matrixd<6, 3> calc_acc(const vector<Transf> & tf0_n_vec,
                             const vector<Matrix3d> & jac_zeta) const;

    vector<Vector2d> _deltaQ;
    Transf _zetaPrior;
    Matrix6d _A;

};

