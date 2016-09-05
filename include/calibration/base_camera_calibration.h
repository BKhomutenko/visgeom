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

template<template<typename> class Projector>
class ExtrinsicCameraCalibration : GenericCameraCalibration<Projector>
{
private:
    vector<double> intrinsicVec;
    vector<CalibrationData> monoCalibDataVec;
    
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

    bool initialize(const string & infoFileName)
    {
        cout << "### Initialize calibration data ###" << endl;
        if (not initializeExtrinsic(infoFileNameStereo)) return false;
        constructGrid();
        return true;
    } 
    
    bool initializeExtrinsic(const string & infoFileName)
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
        
        int inrinsicSize;
        paramFile >> inrinsicSize;
        intrinsicVec.resize(intrinsicSize);
        cout << "Camera model parameters :" << endl;
        for (auto & p: intrinsicVec) 
        {
            paramFile >> p;
            cout << setw(10) << p;
        }
        cout << endl;
        paramFile.ignore();
    
        stereoCalibDataVec.clear();
        
        string imageFolder;
        string imageInfo;    
        string imageName;
        array<double, 6> robotPose;
        
        getline(calibInfoFile, imageFolder);
        while (getline(calibInfoFile, imageName))
        {
            istringstream imageStream(imageInfo);
            
            imageStream >> imageName;
            for (auto & x : robotPose) imageStream >> x;
            
            CalibrationData calibData;
            Vector2dVec point1Vec, point2Vec;
            
            calibData.fileName1 = imageFolder + imageName;
            if (not extractGridProjection(calibData.fileName1,
                    calibData.projection1, checkExtraction)) continue;
            
            calibDataVec.push_back(calibData);
            
            cout << "." << flush;
        }
        cout << "done" << endl;
        return true;
    }
    

    void addStereoResidual(Problem & problem, vector<double> & intrinsic, 
            array<double, 6> & extrinsic,
            const Vector2dVec & projection,
            array<double, 6> & gridExtrinsic)
    {
        typedef DynamicAutoDiffCostFunction<StereoGridProjection<Projector>> stereoProjectionCF;
        StereoGridProjection<Projector> * stereoProjection;
        stereoProjection = new StereoGridProjection<Projector>(projection, grid);
        stereoProjectionCF * costFunction = new stereoProjectionCF(stereoProjection);
        costFunction->AddParameterBlock(intrinsic.size());
        costFunction->AddParameterBlock(6);
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), intrinsic.data(),
                gridExtrinsic.data(), extrinsic.data());
    }
    
    void estimateInitialExtrinsic(const vector<double> & intrinsic, array<double, 6> & extrinsic,
            vector<StereoCalibrationData> & calibDataVec)
    {
        Problem problem;
        for (int i = 0; i < calibDataVec.size(); i++)
        {
            
            typedef DynamicAutoDiffCostFunction<StereoEstimate<Projector>> dynamicProjectionCF;

            StereoEstimate<Projector> * stereoEstimate;
            stereoEstimate = new StereoEstimate<Projector>(calibDataVec[i].projection2, grid,
                                                intrinsic, calibDataVec[i].extrinsic);
            dynamicProjectionCF * costFunction = new dynamicProjectionCF(stereoEstimate);
            costFunction->AddParameterBlock(6);
            costFunction->SetNumResiduals(2 * Nx * Ny);
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), extrinsic.data());   
            
            //run the solver
            
        }
        Solver::Options options;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
    }

    bool estimateExtrinsic(const vector<double> & intrinsic1, 
            const vector<double> & intrinsic2,
            array<double, 6> & extrinsic)
    {
        for (auto & item : monoCalibDataVec1)
        {
            estimateInitialGrid(intrinsic1, item.projection, item.extrinsic);
        }
        for (auto & item : monoCalibDataVec2)
        {
            estimateInitialGrid(intrinsic1, item.projection, item.extrinsic);
        }
        for (auto & item : stereoCalibDataVec)
        {
            estimateInitialGrid(intrinsic1, item.projection1, item.extrinsic);
        }
        estimateInitialExtrinsic(intrinsic2, extrinsic, stereoCalibDataVec);
    }
    
    bool compute(vector<double> & intrinsic1,
            vector<double> & intrinsic2,
            Transformation<double> & transfo)
    {

        array<double, 6> extrinsic = transfo.toArray();
        
        
        cout << "### Extrinsic parameters calibration ###" << endl;
        if (monoCalibDataVec1.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
           
        if (monoCalibDataVec2.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
        
        if (stereoCalibDataVec.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
        
        cout << "Initial board poses estimation " << endl;
        estimateExtrinsic(intrinsic1, intrinsic1, extrinsic);
        
        cout << "Global optimization" << endl;
        Problem problem;
        // Intrinsic init
        for (auto & item : monoCalibDataVec1)
        {
            addIntrinsicResidual(problem, intrinsic1, item.projection, item.extrinsic);
        }
        for (auto & item : monoCalibDataVec2)
        {
            addIntrinsicResidual(problem, intrinsic2, item.projection, item.extrinsic);
        }
            
        // Extrinsic init
        for (auto & item : stereoCalibDataVec)
        {
            addIntrinsicResidual(problem, intrinsic1, item.projection1, item.extrinsic);
            addStereoResidual(problem, intrinsic2, extrinsic, item.projection2, item.extrinsic);
        }
        
        Solver::Options options;
        options.max_num_iterations = 500;
        
        options.function_tolerance = 1e-15;
        options.gradient_tolerance = 1e-10;
        options.parameter_tolerance = 1e-10;
//        options.minimizer_progress_to_stdout = true;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
//        cout << summary.BriefReport() << endl;
        
        transfo = Transformation<double>(extrinsic.data());
    } 
    
    bool residualAnalysis(const vector<double> & intrinsic1, 
            const vector<double> & intrinsic2,
            const Transformation<double> & TleftRight)
    {
        residualAnalysis(intrinsic1, monoCalibDataVec1);
        residualAnalysis(intrinsic2, monoCalibDataVec2);
        vector<CalibrationData> stereoLeftData, stereoRightData;
        for (auto & x : stereoCalibDataVec)
        {
            stereoLeftData.push_back(x.getLeftMono());
            stereoRightData.push_back(x.getRightMono(TleftRight));
        }
        residualAnalysis(intrinsic1, stereoLeftData);
        residualAnalysis(intrinsic2, stereoRightData);
    }
};

