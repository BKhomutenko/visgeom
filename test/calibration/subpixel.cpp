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

#include "io.h"
#include "ocv.h"
#include "eigen.h"

#include "calibration/corner_detector.h"

int main(int argc, char** argv) {

    const int Nx = 8;
    const int Ny = 5;
    
    vector<cv::Point2f> centers;
    Size patternSize(Nx, Ny);
    Vector2dVec cornerVec;
    Mat8u frame = imread(argv[1], 0);
    
//    resize(frame, frame, Size(0, 0), 2, 2);
    
    Mat8u corners1, corners2;
    frame.copyTo(corners1);
    frame.copyTo(corners2);
    bool patternIsFound = findChessboardCorners(frame, patternSize, centers, CV_CALIB_CB_ADAPTIVE_THRESH);
    if (not patternIsFound)
    {
        cout << argv[1] << " : ERROR, pattern not found" << endl;
        return 0;
    }

    
    
        
    cornerVec.resize(Nx * Ny);
    for (int i = 0; i < Nx * Ny; i++)
    {
        cornerVec[i] = Vector2d(centers[i].x, centers[i].y);
    }
    
    double minDist = findMinDistance(cornerVec, Ny, Nx);
    
    CornerDetector detector(frame, min(minDist / 3., 10.));
    
    cout << "before : " << endl;
    for (auto & pt : cornerVec)
    {
        cout << pt.transpose() << endl;
    }
    
    
    
    detector.improveCorners(cornerVec);
    
    cout << "after : " << endl;
    for (auto & pt : cornerVec)
    {
        cout << pt.transpose() << endl;
    }
    
    drawChessboardCorners(corners1, patternSize, Mat(centers), patternIsFound);
    imshow("corners1", corners1);
    
    vector<cv::Point2f> centers2;
    for (auto & x : cornerVec)
    {
        centers2.emplace_back(x[0], x[1]);
    }
    
    drawChessboardCorners(corners2, patternSize, Mat(centers2), patternIsFound);
    imshow("corners2", corners2);
    
    waitKey();
    return 0;
}
