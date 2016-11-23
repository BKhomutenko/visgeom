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
#include "ocv.h"
#include "eigen.h"
#include "ceres.h"

//TODO make a separate header
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <glog/logging.h>

#include "cost_functors.h"

using boost::property_tree::ptree;
using boost::property_tree::read_json;

struct CalibrationData
{
    Vector2dVec projection;
    array<double, 6> extrinsic;
    string fileName;
};

struct CalibrationBoard
{
    string type;
    int nx, ny;
    double size;
};

//TODO choose the model in the parameter file, not as atemplate parameter
template<template<typename> class Projector>
class GenericCameraCalibration
{
protected:
    map<string, double*> intrinsicMap;
    map<string, map<string, double*>> sequenceTransformMap;
    map<string, double*> transformMap;
    Problem globalProblem;
    
    //temporary variables
    int Nx, Ny;
    
public:
    //TODO implement access to the results
    
    ~GenericCameraCalibration()
    {
        for (auto & x : intrinsicMap)
        {
            delete[] x.second;
        }
        for (auto & sequenceMap : sequenceTransformMap)
        {
            for (auto & x : sequenceMap.second)
            {
                delete[] x.second;
            }
        }
        for (auto & x : transformMap)
        {
            delete[] x.second;
        }
    }
    
    bool compute()
    {
               
        //run the solver
        Solver::Options options;
        options.max_num_iterations = 500;
        options.function_tolerance = 1e-10;
        options.gradient_tolerance = 1e-10;
        options.parameter_tolerance = 1e-10;
//        options.minimizer_progress_to_stdout = true;
        Solver::Summary summary;
        Solve(options, &globalProblem, &summary);
        cout << summary.FullReport() << endl;
        
        for (auto & x : intrinsicMap)
        {
            cout << x.first << " : ";
            for (int i = 0; i < Projector<double>::INTRINSIC_COUNT; i++)
            {
                cout << x.second[i] << "  ";
            }
            cout << endl;
        }
    } 
    
    bool addResiduals(const string & infoFileName)
    {
        ptree root;
        read_json(infoFileName, root);
        
        vector<string> transNameVec;
        vector<TransformationStatus> transStatusVec;
        vector<Transformation<double>> constTransVec;
        vector<double> intrinsicVec;
//        vector<double> extrinsicDefault;
        string cameraId, sequenceId;
        
        //read out the transformations
        for (auto & transInfo : root.get_child("transform_chain"))
        {
            transNameVec.emplace_back(transInfo.second.get<string>("name"));
            bool isVariable = transInfo.second.get<bool>("variable");
            bool isDirect = transInfo.second.get<bool>("direct");
            
            //read out the initial value
            vector<double> transform;
            for (auto & x : transInfo.second.get_child("value"))
            {
                transform.push_back( double(x.second.get_value<double>()) );
            }
            
            if (isVariable)
            {
                if (isDirect) transStatusVec.push_back(TRANFORM_DIRECT);
                else transStatusVec.push_back(TRANSFORM_INVERSE);
            }
            else
            {
                //read out the value
                //fill up the const transformation vector
                transStatusVec.push_back(TRANSFORM_CONSTANT);
                constTransVec.emplace_back(transform.data());
                if (not isDirect)
                {
                    constTransVec.back() = constTransVec.back().inverse();
                }
            }
            
            //initialize the transform in the map
            if (transNameVec.back() == "xiImage")
            {
                sequenceId = transInfo.second.get<string>("id");
                if (sequenceTransformMap.find(sequenceId) == sequenceTransformMap.end())
                {
                    sequenceTransformMap[sequenceId] = map<string, double*>();
//                    extrinsicDefault = transform;
                }
            }
            else if (isVariable)
            {
                if (transformMap.find(transNameVec.back()) == transformMap.end())
                {
                    transformMap[transNameVec.back()] = new double[6];
                }
            }
        }
        
        //read out the camera
        cameraId = root.get<string>("camera.name");
        bool isVariableIntrinsic = root.get<bool>("camera.variable");
        if (isVariableIntrinsic)
        {
            if (intrinsicMap.find(sequenceId) == intrinsicMap.end())
            {
                intrinsicMap[cameraId] = new double[Projector<double>::INTRINSIC_COUNT]; //FIXME
                int i = 0;
                for (auto & x : root.get_child("camera.intrinsic"))
                {
                    intrinsicMap[cameraId][i] = double(x.second.get_value<double>());
                    i++;
                }
            }
        }
        else
        {
            for (auto & x : root.get_child("camera.intrinsic"))
            {
                intrinsicVec.push_back( double(x.second.get_value<double>()) );
            }
        }
        //initialize the object
        Nx = root.get<int>("object.cols");
        Ny = root.get<int>("object.rows");
        double sqSize= root.get<double>("object.size");
        Vector3dVec grid;
        grid.reserve(Nx*Ny);
        for (int i = 0; i < Ny; i++)
        {
            for (int j = 0; j < Nx; j++)
            {
               grid.emplace_back(sqSize * j, sqSize * i, 0); 
            }
        }
        
        //initialize the sequence extrinsic, grid projections
        //and add the residuals to the problem
        Vector2dVec projection;
        string prefix = root.get<string>("camera.images.prefix");
        for (auto & x : root.get_child("camera.images.names"))
        {
            //TODO make possible const xiImage
            string filename = x.second.get_value<string>();
            if (not extractGridProjection(prefix + filename, projection,
                    root.get<bool>("params.check_extraction"))) continue;
            
            if (sequenceTransformMap[sequenceId].find(filename) == sequenceTransformMap[sequenceId].end())
            {
                sequenceTransformMap[sequenceId][filename] = new double[6];
                estimateInitialGrid(intrinsicMap[cameraId], projection, grid,
                    sequenceTransformMap[sequenceId][filename]);
            }
            
            // make the vector of pointers
            vector<double*> ptrVec;
            for (int i = 0; i < transNameVec.size(); i++)
            {
                const string & name = transNameVec[i];
                if (name == "xiImage")
                {
                    ptrVec.push_back(sequenceTransformMap[sequenceId][filename]);
                }
                else if (transStatusVec[i] != TRANSFORM_CONSTANT)
                {
                    ptrVec.push_back(transformMap[name]);
                }
            }
            
            
            //add a residual
            typedef DynamicAutoDiffCostFunction<GenericProjection<Projector>> projectionCF;
            GenericProjection<Projector> * boardProjection;
            boardProjection = new GenericProjection<Projector>(projection, grid,
                                        transStatusVec, constTransVec, intrinsicVec);
            projectionCF * costFunction = new projectionCF(boardProjection);
            costFunction->SetNumResiduals(2 * Nx * Ny);
            //add parameter blocks for all nonconstant transformations
            for (int i = 0; i < ptrVec.size(); i++)
            {
                costFunction->AddParameterBlock(6); 
            }
            
            if (isVariableIntrinsic)
            {
                costFunction->AddParameterBlock(Projector<double>::INTRINSIC_COUNT);
                switch (ptrVec.size())
                {
                case 1:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], intrinsicMap[cameraId]);
                    break;
                case 2:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], intrinsicMap[cameraId]);
                    break;
                case 3:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2], intrinsicMap[cameraId]);
                    break;
                case 4:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3], intrinsicMap[cameraId]);
                    break;
                case 5:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3], ptrVec[4], intrinsicMap[cameraId]);
                    break;
                }
            }
            else
            {
                switch (ptrVec.size())
                {
                case 1:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0]);
                    break;
                case 2:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1]);
                    break;
                case 3:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2]);
                    break;
                case 4:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3]);
                    break;
                case 5:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3], ptrVec[4]);
                    break;
                }
            }
        }
        
    }

    bool extractGridProjection(const string & fileName, Vector2dVec & projection, bool checkExtraction)
    {
        Size patternSize(Nx, Ny);
        Mat frame = imread(fileName, 0);
        if (frame.empty())
        {
            cout << fileName << " : ERROR, file not found" << endl;
            return false;
        }
        vector<cv::Point2f> centers;
        bool patternIsFound = findChessboardCorners(frame, patternSize, centers, CV_CALIB_CB_ADAPTIVE_THRESH);
        if (not patternIsFound)
        {
            cout << fileName << " : ERROR, pattern not found" << endl;
            return false;
        }

        if (checkExtraction)
        {
            drawChessboardCorners(frame, patternSize, Mat(centers), patternIsFound);
            imshow("corners", frame);
            char key = waitKey();
            if (key == 'n' or key == 'N')
            {
                cout << fileName << " : ERROR, pattern not accepted" << endl;
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

    void estimateInitialGrid(const double * const intrinsic,
            const Vector2dVec & projection, const Vector3dVec & grid,
            double * const extrinsic)
    {
        Problem problem;
        typedef DynamicAutoDiffCostFunction<GridEstimate<Projector>> dynamicProjectionCF;

        GridEstimate<Projector> * boardEstimate;
        
        fill(extrinsic, extrinsic + 6, 0);
        extrinsic[2] = 1;
        
        //compute initial orientation
        auto v = projection[1] - projection[0];
        float alpha = atan2(v[1], v[0]);
        extrinsic[5] = alpha;
        
        //optimize the position
        vector<double> intrinsicVec(intrinsic, intrinsic + Projector<double>::INTRINSIC_COUNT);
        boardEstimate = new GridEstimate<Projector>(projection,
                                    grid, intrinsicVec);
        dynamicProjectionCF * costFunction = new dynamicProjectionCF(boardEstimate);
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), extrinsic);

        //run the solver
        Solver::Options options;
        options.max_num_iterations = 250;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
    }
    
/*


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
    }*/
};

