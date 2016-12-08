/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUdouble ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 


#pragma once

#include "eigen.h"

#include "geometry/geometry.h"

enum CameraIdx {CAMERA_1, CAMERA_2};

const double DEFAULT_DEPTH = 5;
const double MIN_DEPTH = 0.1;
const double MAX_DEPTH = 100;
const double DEFAULT_SIGMA_DEPTH = 30;
const double DEFAULT_COST_DEPTH = 5;
const double COST_CHANGE = 1;
const double OUT_OF_RANGE = 0.0;
const double SIGMA_COEFF = 1 / 1.7;

/*

triangulation is done by solving

{ p (l1 p - l2 q - t) = 0
{ q (l1 p - l2 q - t) = 0

p, q are direction vectors;
t is a translation taken from the transformation

jac is a kinematic jacobian  l_dot = jac  *   (  v  )
                                              (omega)
Memory must be allocated
*/
void triangulate(const Transf xi, 
        const Vector3dVec & xVec1, const Vector3dVec & xVec2,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL);

// num / denom if denom is big enough
// linear extrapolation otherwise
double regDiv(double num, double denom, double EPS);

/*
Triangulation is done by solving:

{ t.(l1 p - l2 q - t) = 0
{ (p + q).(l1 p - l2 q - t) = 0

p, q are direction vectors;
t is the translation taken from the transformation

in case if p and q are close to parallel, the reconstruction is regularized
*/
void triangulateRegular(const Transf xi, 
        const Vector3dVec & xVec1, const Vector3dVec & xVec2,
        double * res1, double * res2 = NULL,
        double * jac1 = NULL, double * jac2 = NULL);

/*
dynamic algorithm for descriptor comparison

returns a verctor of costs for each possible disparity value
*/
vector<int> compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec, int flawCost);
        
        
