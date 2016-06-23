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

#include "generic_calibration.h"

struct RobotCalibrationData
{
    Vector2dVec projection;
    array<double, 6> robotPose;
    string fileName;
};

template<template<typename> class Projector>
struct TransformedGridProjection 
{
    TransformedGridProjection(const Vector2dVec & proj,
            const Vector3dVec & grid, 
            const Transformation<double> & TterminalGrid)
    : _proj(proj), _grid(grid), _TterminalGrid(TterminalGrid) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> TcameraTerminal(params[1]);
        Transformation<T> TterminalGrid = _TterminalGrid.template cast<T>();
        Vector3Vec<T> transformedPoints(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            transformedPoints[i] = _grid[i].template cast<T>();
        }
        TcameraTerminal.compose(TterminalGrid).transform(transformedPoints, transformedPoints);
        
        Projector<T> projector;
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(params[0], transformedPoints[i].data(), modProj.data())) 
            {
                Vector2<T> diff = _proj[i].template cast<T>() - modProj;
                residual[2*i] = diff[0];
                residual[2*i + 1] = diff[1];
            }
            else
            {
                residual[2*i] = T(0.);
                residual[2*i + 1] = T(0.);
            }
        }
        return true;
    }
    
    const Transformation<double> _TterminalGrid;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
};

template<template<typename> class Projector>
class EndEffectorCameraCalibration : GenericCameraCalibration<Projector>
{
private:
    vector<RobotCalibrationData> robotCalibDataVec;
    array<double, 6> gridPose;
    
    // methods
    using GenericCameraCalibration<Projector>::initializeIntrinsic;
    using GenericCameraCalibration<Projector>::extractGridProjection;
    using GenericCameraCalibration<Projector>::constructGrid;
    using GenericCameraCalibration<Projector>::estimateInitialGrid;
    using GenericCameraCalibration<Projector>::addIntrinsicResidual;
    using GenericCameraCalibration<Projector>::residualAnalysis;
    
    // fields
    using GenericCameraCalibration<Projector>::grid;
    using GenericCameraCalibration<Projector>::Nx;
    using GenericCameraCalibration<Projector>::Ny;
    using GenericCameraCalibration<Projector>::sqSize;
    using GenericCameraCalibration<Projector>::outlierThresh;
        
public:


    //TODO chanche the file formatting
    bool initialize(const string & infoFileName)
    {
        cout << "### Initialize calibration data ###" << endl;
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
        
        for (auto & x : gridPose) calibInfoFile >> x;  // grid in robot's base frame
        calibInfoFile.ignore();
        
        robotCalibDataVec.clear();
        string imageFolder;
        string imageInfo;
        getline(calibInfoFile, imageFolder);
        while (getline(calibInfoFile, imageInfo))
        {
            RobotCalibrationData calibData;
            istringstream imageStream(imageInfo);
            
            string imageName;
            imageStream >> imageName;
            for (auto & x : calibData.robotPose) imageStream >> x;
                        
            calibData.fileName = imageFolder + imageName;
            bool isExtracted = extractGridProjection(calibData.fileName,
                    calibData.projection, checkExtraction);

            if (not isExtracted)
            {
                continue;
            }

            robotCalibDataVec.push_back(calibData);

            cout << "." << flush;
        }
        cout << "done" << endl;
        constructGrid();
        return true;
    }
    
    bool compute(vector<double> & intrinsic, Transformation<double> & TterminalCamera)
    {
        array<double, 6> endEffectorPose;
        cout << "### Initial extrinsic estimation ###" << endl;

        if (robotCalibDataVec.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
        
        //initial end effector transform estimation
        estimateEndEffectorTransform(intrinsic, endEffectorPose, robotCalibDataVec[0]); 

        cout << "### Intrinsic parameters calibration ###" << endl;
        // Problem initialization
        Problem problem;
        initEndEffectorProblem(problem, intrinsic, endEffectorPose, robotCalibDataVec);
        //run the solver
        Solver::Options options;
        options.max_num_iterations = 500;
        options.function_tolerance = 1e-10;
        options.gradient_tolerance = 1e-10;
        options.parameter_tolerance = 1e-10;
//        options.minimizer_progress_to_stdout = true;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
        cout << summary.FullReport() << endl;
        TterminalCamera = Transformation<double>(endEffectorPose.data()).inverse();
    } 
    
    void estimateEndEffectorTransform(const vector<double> & intrinsic, 
            array<double, 6> & endEffectorPose,
            const RobotCalibrationData & item)
    {
        Problem problem;
        typedef DynamicAutoDiffCostFunction<GridEstimate<Projector>> dynamicProjectionCF;

        GridEstimate<Projector> * boardEstimate;
        
        //compute initial orientation
        array<double, 6> boardPose = {0,  0, 1, 0,  0, 0};
        auto v = item.projection[2] - item.projection[0];
        float alpha = atan2(v[1], v[0]);
        boardPose[5] = alpha;
        
        //optimize the position
        boardEstimate = new GridEstimate<Projector>(item.projection,
                                    grid, intrinsic);
        dynamicProjectionCF * costFunction = new dynamicProjectionCF(boardEstimate);
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), boardPose.data());

        //run the solver
        Solver::Options options;
        options.max_num_iterations = 50;
        Solver::Summary summary;
        Solve(options, &problem, &summary);

        // rough estimation of the end effector transformation
        Transformation<double> TbaseGrid(gridPose.data());
        Transformation<double> TbaseTerminal(item.robotPose.data());
        Transformation<double> TgridTerminal = TbaseGrid.inverseCompose(TbaseTerminal);
        Transformation<double> TcameraGrid(boardPose.data());
        Transformation<double> TcameraTerminal = TcameraGrid.compose(TgridTerminal);
        endEffectorPose = TcameraTerminal.toArray();
    }
    
    void initEndEffectorProblem(Problem & problem, vector<double> & intrinsic,
            array<double, 6> & endEffectorPose,
            vector<RobotCalibrationData> & calibDataVec)
    {
        typedef DynamicAutoDiffCostFunction<TransformedGridProjection<Projector>> projectionCF;
        Transformation<double> TbaseGrid(gridPose.data());
        for (auto & item : calibDataVec)
        {
            Transformation<double> TbaseTerminal(item.robotPose.data());
            Transformation<double> TterminalGrid = TbaseTerminal.inverseCompose(TbaseGrid);
            
            TransformedGridProjection<Projector> * gridProjection =
                    new TransformedGridProjection<Projector>(item.projection, grid, TterminalGrid);
            projectionCF * costFunction = new projectionCF(gridProjection);
            costFunction->AddParameterBlock(intrinsic.size());
            costFunction->AddParameterBlock(6);
            costFunction->SetNumResiduals(2 * Nx * Ny);
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), intrinsic.data(),
                    endEffectorPose.data());
        }
    }
    
    void residualAnalysis(const vector<double> & intrinsic, Transformation<double> & TterminalCamera)
    {
        vector<CalibrationData> monoCalibDataVec;
        Transformation<double> TbaseGrid(gridPose.data());
        for (auto & item : robotCalibDataVec)
        {
            Transformation<double> TbaseTerminal(item.robotPose.data());
            Transformation<double> TterminalGrid = TbaseTerminal.inverseCompose(TbaseGrid);
            CalibrationData calibData;
            calibData.fileName = item.fileName;
            calibData.extrinsic = TterminalCamera.inverseCompose(TterminalGrid).toArray();
            calibData.projection = item.projection;
            monoCalibDataVec.push_back(calibData);
        }
        residualAnalysis(intrinsic, monoCalibDataVec);
    }

};

