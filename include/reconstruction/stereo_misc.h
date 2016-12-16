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

enum CameraIdx {CAMERA_1, CAMERA_2};

const double DEFAULT_DEPTH = 5;
const double MIN_DEPTH = 0.25;
const double MAX_DEPTH = 100;
const double DEFAULT_SIGMA_DEPTH = 30;
const double DEFAULT_COST_DEPTH = 5;
const double COST_CHANGE = 1;
const double OUT_OF_RANGE = 0.0;
const double SIGMA_COEFF = 1 / 1.7;


