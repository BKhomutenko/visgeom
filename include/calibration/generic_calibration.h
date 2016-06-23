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

#include <cmath>
#include <algorithm>

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "cost_functors.h"

using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::SoftLOneLoss;

int nObjects = 0;

struct CalibrationData
{
    Vector2dVec projection;
    array<double, 6> extrinsic;
    string fileName;
};

template<template<typename> class Projector>
class GenericCameraCalibration
{
protected:
    int Nx, Ny;
    double sqSize;
    double outlierThresh;
    Vector3dVec grid;
public:

    //TODO chanche the file formatting
    bool initializeIntrinsic(const string & infoFileName,
            vector<CalibrationData> & calibDataVec)
    {
        // open the file and read the data
        ifstream calibInfoFile(infoFileName);
        if (not calibInfoFile.is_open())
        {
            cout << infoFileName << " : ERROR, file is not found" << endl;
            return false;
        }
        bool checkExtraction;
        calibInfoFile >> Nx >> Ny >> sqSize >> outlierThresh >> checkExtraction;
        calibInfoFile.ignore();  // To get to the next line

        calibDataVec.clear();
        string imageFolder;
        string imageName;
        getline(calibInfoFile, imageFolder);
        while (getline(calibInfoFile, imageName))
        {
            CalibrationData calibData;
            bool isExtracted;

            calibData.fileName = imageFolder + imageName;
            isExtracted = extractGridProjection(calibData.fileName, calibData.projection, checkExtraction);

            if (not isExtracted)
            {
                continue;
            }

            calibData.extrinsic = array<double, 6>{0, 0, 1, 0, 0, 0};
            calibDataVec.push_back(calibData);

            cout << "." << flush;
        }
        cout << "done" << endl;
        return true;
    }

    bool extractGridProjection(const string & fileName, Vector2dVec & projection, bool checkExtraction)
    {
        Size patternSize(Nx, Ny);
        Mat frame = imread(fileName, 0);

        vector<cv::Point2f> centers;
        bool patternIsFound = findChessboardCorners(frame, patternSize, centers);
        if (not patternIsFound)
        {
            cout << fileName << " : ERROR, pattern is not found" << endl;
            return false;
        }

        if (checkExtraction)
        {
            drawChessboardCorners(frame, patternSize, Mat(centers), patternIsFound);
            imshow("corners", frame);
            char key = waitKey();
            if (key == 'n' or key == 'N')
            {
                cout << fileName << " : ERROR, pattern is not accepted" << endl;
                return false;
            }
        }

        projection.resize(Nx * Ny);
        for (int i = 0; i < Nx * Ny; i++)
        {
            projection[i] = Vector2d(centers[i].x, centers[i].y);
        }
        return true;
    }

    void constructGrid()
    {
        grid.resize(Nx * Ny);
        for (int i = 0; i < Nx * Ny; i++)
        {
            grid[i] = Vector3d(sqSize * (i % Nx), sqSize * (i / Nx), 0);
        }
    }

    void estimateInitialGrid(const vector<double> & intrinsic,
            const Vector2dVec & projection,
            array<double, 6> & extrinsic)
    {
        Problem problem;
        typedef DynamicAutoDiffCostFunction<GridEstimate<Projector>> dynamicProjectionCF;

        GridEstimate<Projector> * boardEstimate;
        
        //compute initial orientation
        
        auto v = projection[1] - projection[0];
        float alpha = atan2(v[1], v[0]);
        extrinsic[5] = alpha;
        
        //optimize the position
        boardEstimate = new GridEstimate<Projector>(projection,
                                    grid, intrinsic);
        dynamicProjectionCF * costFunction = new dynamicProjectionCF(boardEstimate);
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), extrinsic.data());

        //run the solver
        Solver::Options options;
        options.max_num_iterations = 250;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
    }

    void addIntrinsicResidual(Problem & problem, vector<double> & intrinsic,
            const Vector2dVec & projection,
            array<double, 6> & extrinsic)
    {
        typedef DynamicAutoDiffCostFunction<GridProjection<Projector>> projectionCF;
        GridProjection<Projector> * boardProjection;
        boardProjection = new GridProjection<Projector>(projection, grid);
        projectionCF * costFunction = new projectionCF(boardProjection);
        costFunction->AddParameterBlock(intrinsic.size());
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), intrinsic.data(), extrinsic.data());
    }

    void residualAnalysis(const vector<double> & intrinsic,
            const vector<CalibrationData> & calibDataVec)
    {
        double Ex = 0, Ey = 0;
        double Emax = 0;
        Mat_<float> errorPlot(400, 400, 0.f);
        circle(errorPlot, Point(200, 200), 50, 0.4, 1);
        circle(errorPlot, Point(200, 200), 100, 0.4, 1);
        circle(errorPlot, Point(200, 200), 150, 0.4, 1);
        for (int ptIdx = 0; ptIdx < calibDataVec.size(); ptIdx++)
        {
                Vector3dVec transfModelVec;
                Transformation<double> TcamGrid(calibDataVec[ptIdx].extrinsic.data());
                TcamGrid.transform(grid, transfModelVec);

                Vector2dVec projModelVec(transfModelVec.size());
                Projector<double> projector;
                for (int i = 0; i < transfModelVec.size(); i++)
                {
                    projector(intrinsic.data(),
                            transfModelVec[i].data(),
                            projModelVec[i].data());
                }

                Mat frame = imread(calibDataVec[ptIdx].fileName, 0);

                bool outlierDetected = false;
                for (int i = 0; i < Nx * Ny; i++)
                {
                    Vector2d p = calibDataVec[ptIdx].projection[i];
                    Vector2d pModel = projModelVec[i];
                    Vector2d delta = p - pModel;
                    int y0 = floor(delta(1)*100+ 200), x0 = floor(delta(0)*100+ 200);
                    for (auto y : {y0, y0 + 1})
                    {
                        for (auto x : {x0, x0 + 1})
                        {
                            if (y >= 0 and y < 400 and x >= 0 and x < 400)
                                errorPlot(y, x) += 0.2;
                        }
                    }
                    double dx = delta[0] * delta[0];
                    double dy = delta[1] * delta[1];
                    if (outlierThresh != 0 and dx + dy > outlierThresh * outlierThresh)
                    {
                        outlierDetected = true;
                        cout << calibDataVec[ptIdx].fileName << " # " << i << endl;
                        cout << delta.transpose() << endl;
                        circle(frame, Point(pModel(0), pModel(1)), 8, 105, 3);
                    }
                    else
                    {
                        circle(frame, Point(pModel(0), pModel(1)), 8, 190, 3);
                    }
                    if (dx + dy > Emax)
                    {
                        Emax = dx + dy;
                    }
                    Ex += dx;
                    Ey += dy;
                }
                if (outlierDetected)
                {
                    imshow("reprojection", frame);
                    waitKey();
                }
        }
        Ex /= calibDataVec.size() * Nx * Ny;
        Ey /= calibDataVec.size() * Nx * Ny;
        Ex = sqrt(Ex);
        Ey = sqrt(Ey);
        Emax = sqrt(Emax);
        cout << "Ex = " << Ex << "; Ey = " << Ey << "; Emax = " << Emax << endl;
        imshow("errorPlot", errorPlot);
        waitKey();
    }
};

