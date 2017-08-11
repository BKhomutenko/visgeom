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
#include <opencv2/features2d/features2d.hpp>

// Data types
using cv::Point;
using cv::Scalar;
using cv::Point2d;
using cv::Size;
using cv::Rect;
using cv::Scalar;

using cv::Mat;
using cv::Mat_;
using cv::MatND;
using Mat32f = Mat_<float>;
using Mat8u = cv::Mat_<uint8_t>;
using Mat8uc3 = cv::Mat_<cv::Vec3b>;
using Mat16s = cv::Mat_<int16_t>;
using Mat32s = cv::Mat_<int32_t>;

// Functions
using cv::swap;
using cv::imshow;
using cv::imread;
using cv::imwrite;
using cv::waitKey;
using cv::circle;

//point features
using cv::BRISK;
using cv::KeyPoint;
using cv::DMatch;
using cv::BFMatcher;

//interface
using cv::setMouseCallback;

//TODO put it elsewhere
template<typename T, typename M, typename Q>
T bilinear(const M & src, Q x, Q y)
{
    int u(x), v(y);
    Q dx = x - u;
    Q dy = y - v;
    Q dx2 = 1 - dx;
    
    //check the limits
    bool fail = false;
    if (fail |= u < 0) u = 0;
    else if (fail |= u > src.cols - 2) u = src.cols - 1;
    if (fail |= v < 0) v = 0;
    else if (fail |= v > src.rows - 2) v = src.rows - 1;
    if (fail) return src(v, u);
    
    const T i00 = src(v, u);
    const T i01 = src(v, u + 1);
    const T i10 = src(v + 1, u);
    const T i11 = src(v + 1, u + 1);
    return (i11*dx + i10*dx2)*dy + (i01*dx + i00*dx2)*(1-dy);
}

inline void cross(Mat& img, double x, double y, int size, const Scalar& color,
            int thickness=1, int lineType=8, int shift=0)
{
    line(img, Point(x - size, y - size),
            Point(x + size, y + size),
            color, thickness, lineType, shift);
    line(img, Point(x - size, y + size),
            Point(x + size, y - size),
            color, thickness, lineType, shift);
}

template <typename T>
void drawPoints(Mat& img, T pointVec)
{
    if (pointVec.size() == 0) return;
    if (img.channels() == 1)
    {
        int color = 0;
        const int colorStep = 155 / pointVec.size();
        
        for (auto & pt : pointVec)
        {
            cross(img, pt[0], pt[1], 10, color, 3);
            cross(img, pt[0], pt[1], 9, color + 100, 1);
            color += colorStep;
        }
    }
    else if (img.channels() == 3)
    {
        int color = 255;
        const int colorStep = 255 / pointVec.size();
        
        for (auto & pt : pointVec)
        {
            int color2 = sqrt(255 * 255 - color * color);
            cross(img, pt[0], pt[1], 10, Scalar(color, 0, color2) , 3);
            color -= colorStep;
        }
    }
}

