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

using ceres::SoftLOneLoss;

template<template<typename> class Projector>
struct StereoGridProjection 
{
    StereoGridProjection(const vector<Eigen::Vector2d> & proj,
            const vector<Eigen::Vector3d> & orig) : _proj(proj), _orig(orig) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> TrefGrid(params[1]);
        Transformation<T> TrefCam(params[2]);
        Transformation<T> TcamGrid = TrefCam.inverseCompose(TrefGrid);
        
        vector<Vector3<T>> transformedPoints(_orig.size());
        for (int i = 0; i < _orig.size(); i++)
        {
            transformedPoints[i] = _orig[i].template cast<T>();
        }
        TcamGrid.transform(transformedPoints, transformedPoints);
        
        Projector<T> projector;
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            projector(params[0], transformedPoints[i].data(), modProj.data());
            Vector2<T> diff = _proj[i].template cast<T>() - modProj;
            residual[2*i] = T(diff[0]);
            residual[2*i + 1] = T(diff[1]);
        }
        return true;
    }
    
    const vector<Vector2d> & _proj;
    const vector<Vector3d> & _orig;
};

template<template<typename> class Projector>
struct StereoEstimate
{
    StereoEstimate(const vector<Vector2d> & proj, const vector<Vector3d> & orig,
    const vector<double> & camParams, const array<double, 6> & extrinsic) 
    : _proj(proj), _orig(orig), _camParams(camParams),
      _extrinsic(extrinsic) {}
            
    template <typename T>
    bool operator()(const T * const * params,
                    T* residual) const 
    {
        array<T, 6> extrinsic;
        for (int i = 0; i < _extrinsic.size(); i++)
        {
            extrinsic[i] = T(_extrinsic[i]);
        }
        Transformation<T> TrefGrid(extrinsic.data());
        Transformation<T> TrefCam(params[0]);
        Transformation<T> TcamGrid = TrefCam.inverseCompose(TrefGrid);
        vector<Vector3<T>> transformedPoints(_orig.size());
        for (int i = 0; i < _orig.size(); i++)
        {
            transformedPoints[i] = _orig[i].template cast<T>();
        }
        TcamGrid.transform(transformedPoints, transformedPoints);
        
        vector<T> camParamsT(_camParams.begin(), _camParams.end());
        Projector<T> projector;
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            projector(camParamsT.data(), transformedPoints[i].data(), modProj.data());
            Vector2<T> diff = _proj[i].template cast<T>() - modProj;
            residual[2*i] = T(diff[0]);
            residual[2*i + 1] = T(diff[1]);
        }
        return true;
    }
    
    const array<double, 6> & _extrinsic;
    const vector<double> & _camParams;
    const vector<Vector2d> & _proj;
    const vector<Vector3d> & _orig;
};

template<template<typename> class Projector>
class ExtrinsicCameraCalibration : GenericCameraCalibration<Projector>
{
private:
    vector<CalibrationData> monoCalibDataVec1;
    vector<CalibrationData> monoCalibDataVec2;
    vector<CalibrationData> stereoCalibDataVec1;
    vector<CalibrationData> stereoCalibDataVec2;
    
    // methods
    using GenericCameraCalibration<Projector>::initializeIntrinsic;
    using GenericCameraCalibration<Projector>::extractGridProjection;
    using GenericCameraCalibration<Projector>::constructGrid;
    using GenericCameraCalibration<Projector>::estimateInitialGrid;
    using GenericCameraCalibration<Projector>::initIntrinsicProblem;
    using GenericCameraCalibration<Projector>::residualAnalysis;
    
    // fields
    using GenericCameraCalibration<Projector>::grid;
    using GenericCameraCalibration<Projector>::Nx;
    using GenericCameraCalibration<Projector>::Ny;
    using GenericCameraCalibration<Projector>::sqSize;
    using GenericCameraCalibration<Projector>::outlierThresh;
    
public:

    bool initialize(const string & infoFileName1, const string & infoFileName2, 
            const string & infoFileNameStereo)
    {
        cout << "### Initialize calibration data ###" << endl;
        if (not initializeIntrinsic(infoFileName1, monoCalibDataVec1))
        {
            return false;
        }
        if (not initializeIntrinsic(infoFileName2, monoCalibDataVec2))
        {
            return false;
        }
        constructGrid();
        
        if (not initializeExtrinsic(infoFileNameStereo)) return false;
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
        
        stereoCalibDataVec1.clear();
        stereoCalibDataVec2.clear();
        
        string imageFolder;
        string imageName;    
        string leftPref, rightPref;
        
        getline(calibInfoFile, imageFolder);
        getline(calibInfoFile, leftPref);
        getline(calibInfoFile, rightPref);
        while (getline(calibInfoFile, imageName))
        {
            CalibrationData calibDataLeft, calibDataRight;
            vector<Vector2d> point1Vec, point2Vec;
            bool isExtracted1, isExtracted2;
            
            calibDataLeft.fileName = imageFolder + leftPref + imageName;
            isExtracted1 = extractGridProjection(calibDataLeft, checkExtraction);
            
            calibDataRight.fileName = imageFolder + rightPref + imageName;                                     
            isExtracted2 = extractGridProjection(calibDataRight, checkExtraction);
            
            if (not isExtracted1 or not isExtracted2) 
            {
                continue;
            }
            
            calibDataLeft.extrinsic = ArraySharedPtr(new array<double, 6>{0, 0, 1, 0, 0, 0});
            calibDataRight.extrinsic = calibDataLeft.extrinsic;
            
            stereoCalibDataVec1.push_back(calibDataLeft);
            stereoCalibDataVec2.push_back(calibDataRight);
            
            cout << "." << flush;
        }
        cout << "done" << endl;
        return true;
    }
    

    void initStereoProblem(Problem & problem, vector<double> & intrinsic, 
            array<double, 6> & extrinsic,
            vector<CalibrationData> & calibDataVec)
    {
        typedef DynamicAutoDiffCostFunction<StereoGridProjection<Projector>> stereoProjectionCF;
        for (unsigned int i = 0; i < calibDataVec.size(); i++)
        {
            StereoGridProjection<Projector> * stereoProjection;
            stereoProjection = new StereoGridProjection<Projector>(
                                        calibDataVec[i].projection,
                                        grid);
            
            stereoProjectionCF * costFunction = new stereoProjectionCF(stereoProjection);
            costFunction->AddParameterBlock(intrinsic.size());
            costFunction->AddParameterBlock(6);
            costFunction->AddParameterBlock(6);
            costFunction->SetNumResiduals(2 * Nx * Ny);
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), intrinsic.data(),
                    calibDataVec[i].extrinsic->data(), extrinsic.data());
            
        }
    }
    
    void estimateInitialExtrinsic(const vector<double> & intrinsic, array<double, 6> & extrinsic,
            vector<CalibrationData> & calibDataVec)
    {
        Problem problem;
        for (int i = 0; i < calibDataVec.size(); i++)
        {
            
            typedef DynamicAutoDiffCostFunction<StereoEstimate<Projector>> dynamicProjectionCF;

            StereoEstimate<Projector> * stereoEstimate;
            stereoEstimate = new StereoEstimate<Projector>(calibDataVec[i].projection, grid,
                                                intrinsic, *(calibDataVec[i].extrinsic));
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
        estimateInitialGrid(intrinsic1, monoCalibDataVec1);
        estimateInitialGrid(intrinsic2, monoCalibDataVec2);
        estimateInitialGrid(intrinsic1, stereoCalibDataVec1);
        estimateInitialExtrinsic(intrinsic2, extrinsic, stereoCalibDataVec2);
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
        
        if (stereoCalibDataVec1.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
        
        cout << "Initial board poses estimation " << endl;
        estimateExtrinsic(intrinsic1, intrinsic1, extrinsic);
        
        cout << "Global optimization" << endl;
        Problem problem;
        // Intrinsic init
        initIntrinsicProblem(problem, intrinsic1, monoCalibDataVec1);
        initIntrinsicProblem(problem, intrinsic2, monoCalibDataVec2);
            
        // Extrinsic init
        initIntrinsicProblem(problem, intrinsic1, stereoCalibDataVec1);
        initStereoProblem(problem, intrinsic2, extrinsic, stereoCalibDataVec2);
            
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
            const Transformation<double> & transfo)
    {
        residualAnalysis(intrinsic1, monoCalibDataVec1);
        residualAnalysis(intrinsic2, monoCalibDataVec2);
        residualAnalysis(intrinsic1, stereoCalibDataVec1);
        residualAnalysis(intrinsic2, stereoCalibDataVec2, transfo);
    }
};

