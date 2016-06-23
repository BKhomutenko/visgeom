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
OpenCV include and using
*/

#pragma once

#include <opencv2/opencv.hpp>

// Data types
using cv::Point;
using cv::Point2d;
using cv::Size;
using cv::Rect;
using cv::Scalar;

using cv::Mat;
using cv::Mat_;
using cv::MatND;
using Mat32f = Mat_<float>;
using Mat8u = cv::Mat_<uint8_t>;
using Mat16s = cv::Mat_<int16_t>;
using Mat32s = cv::Mat_<int32_t>;

// Functions
using cv::imshow;
using cv::imread;
using cv::imwrite;
using cv::waitKey;

using cv::circle;
