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

#include <string> 

using ceres::SoftLOneLoss;

template<template<typename> class Projector>
struct StereoGridProjection 
{
    StereoGridProjection(const Vector2dVec & proj,
            const Vector3dVec & orig) : _proj(proj), _orig(orig) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> TrefGrid(params[1]);
        Transformation<T> TrefCam(params[2]);
        Transformation<T> TcamGrid = TrefCam.inverseCompose(TrefGrid);
        
        Vector3Vec<T> transformedPoints(_orig.size());
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
    
    const Vector2dVec _proj;
    const Vector3dVec _orig;
};

template<template<typename> class Projector>
struct StereoEstimate
{
    StereoEstimate(const Vector2dVec & proj, const Vector3dVec & orig,
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
        Vector3Vec<T> transformedPoints(_orig.size());
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
    
    const array<double, 6> _extrinsic;
    const vector<double> _camParams;
    const Vector2dVec _proj;
    const Vector3dVec _orig;
};


struct StereoCalibrationData
{
    Vector2dVec projection1;
    Vector2dVec projection2;
    array<double, 6> extrinsic;
    string fileName1, fileName2;
    
    CalibrationData getLeftMono()
    {
        CalibrationData res;
        res.projection = projection1;
        res.fileName = fileName1;
        res.extrinsic = extrinsic;
        return res;
    }
    
    CalibrationData getRightMono(const Transformation<double> & TleftRight)
    {
        CalibrationData res;
        res.projection = projection2;
        res.fileName = fileName2;
        Transformation<double> TleftGrid(extrinsic.data());
        Transformation<double> TcameraGrid = TleftRight.inverseCompose(TleftGrid);
        res.extrinsic = TcameraGrid.toArray();
        return res;
    }
};

template<template<typename> class Projector>
class ExtrinsicCameraCalibration : GenericCameraCalibration<Projector>
{
private:
    vector<CalibrationData> monoCalibDataVec1;
    vector<CalibrationData> monoCalibDataVec2;
    vector<StereoCalibrationData> stereoCalibDataVec;
    
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

    bool initialize(const string & infoFileName1, const string & infoFileName2, 
            const string & infoFileNameStereo)
    {
        cout << "### Initialize calibration data ###" << endl;
        for (auto & x : {"ololo", "trololo"})
        {
            cout << x << endl;
        }
        if (not initializeIntrinsic(infoFileName1, monoCalibDataVec1))
        {
            return false;
        }
        if (not initializeIntrinsic(infoFileName2, monoCalibDataVec2))
        {
            return false;
        }
        constructGrid();
        
        vector<double>{1, 2, 1.2};
        
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
        
        stereoCalibDataVec.clear();
        
        string imageFolder;
        string imageName;    
        string leftPref, rightPref;
        
        getline(calibInfoFile, imageFolder);
        getline(calibInfoFile, leftPref);
        getline(calibInfoFile, rightPref);
        while (getline(calibInfoFile, imageName))
        {
            StereoCalibrationData stereoCalibData;
            Vector2dVec point1Vec, point2Vec;
            bool isExtracted1, isExtracted2;
            
            stereoCalibData.fileName1 = imageFolder + leftPref + imageName;
            isExtracted1 = extractGridProjection(stereoCalibData.fileName1,
                    stereoCalibData.projection1, checkExtraction);
            
            stereoCalibData.fileName2 = imageFolder + rightPref + imageName;                                     
            isExtracted2 = extractGridProjection(stereoCalibData.fileName2,
                    stereoCalibData.projection2, checkExtraction);
            
            if (not isExtracted1 or not isExtracted2) 
            {
                continue;
            }
            
            stereoCalibData.extrinsic = array<double, 6>{0, 0, 1, 0, 0, 0};
            stereoCalibDataVec.push_back(stereoCalibData);
            
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

