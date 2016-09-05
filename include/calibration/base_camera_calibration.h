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
struct BaseGridProjection 
{
    BaseGridProjection(const Vector2dVec & proj,
            const Vector3dVec & grid, 
            const Transformation<double> & TorigBase,
            const vector<double> & intrinsicVec)
    : _proj(proj), _grid(grid), _TbaseOrig(TorigBase.inverse()), _intrinsicVec(intrinsicVec) {}
            
    template <typename T>
    bool operator()(const T * const* params,
                    T* residual) const 
    {
        Transformation<T> TcameraBase(params[0]);
        Transformation<T> TorigGrid(params[1]);
        Transformation<T> TbaseOrig = _TbaseOrig.template cast<T>();
        Vector3Vec<T> transformedPoints(_grid.size());
        for (int i = 0; i < _grid.size(); i++)
        {
            transformedPoints[i] = _grid[i].template cast<T>();
        }
        Transformation<T> TcameraGrid = TcameraBase.compose(TbaseOrig).compose(TorigGrid);
        TcameraGrid.transform(transformedPoints, transformedPoints);
        
        Projector<T> projector;
        vector<T> intrinsicTVec;
        intrinsicTVec.reserve(_intrinsicVec.size());
        for (auto & x : _intrinsicVec)
        {
            intrinsicTVec.emplace_back(x);
        }
        for (unsigned int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2<T> modProj;
            if (projector(intrinsicTVec.data(), transformedPoints[i].data(), modProj.data())) 
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
    
    const Transformation<double> _TbaseOrig;
    const Vector2dVec _proj;
    const Vector3dVec _grid;
    const vector<double> _intrinsicVec;
};

template<template<typename> class Projector>
class BaseTransformationCalibration : GenericCameraCalibration<Projector>
{
private:
    vector<double> intrinsicVec;
    vector<CalibrationData> baseCalibDataVec;
    array<double, 6> initialBaseCameraGuess; //FIXME make a transformation
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
        if (not initializeExtrinsic(infoFileName)) return false;
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
        
        int intrinsicSize;
        calibInfoFile >> intrinsicSize;
        intrinsicVec.resize(intrinsicSize);
        cout << "Camera model parameters :" << endl;
        for (auto & p: intrinsicVec) 
        {
            calibInfoFile >> p;
            cout << setw(10) << p;
        }
        cout << endl;
        calibInfoFile.ignore();
    
        cout << "Initial extrinsic guess :" << endl;
        for (auto & p: initialBaseCameraGuess) 
        {
            calibInfoFile >> p;
            cout << setw(10) << p;
        }
        cout << endl;
        calibInfoFile.ignore();
        
        baseCalibDataVec.clear();
        
        string imageFolder;
        string imageInfo;    
        string imageName;
        array<double, 6> robotPose;
        
        getline(calibInfoFile, imageFolder);
        while (getline(calibInfoFile, imageInfo))
        {
            istringstream imageStream(imageInfo);
            
            CalibrationData calibData;
            
            imageStream >> imageName;
            for (auto & x : calibData.extrinsic) imageStream >> x;
            
            calibData.fileName = imageFolder + imageName;
            if (not extractGridProjection(calibData.fileName,
                    calibData.projection, checkExtraction)) continue;
            
            baseCalibDataVec.push_back(calibData);
            
            cout << "." << flush;
        }
        cout << "done" << endl;
        return true;
    }
    

    void addBaseResidual(Problem & problem, const vector<double> & intrinsic, 
            const array<double, 6> & poseOrigBase, const Vector2dVec & projection,
            array<double, 6> & poseCameraBase, array<double, 6> & poseOrigGrid)
    {
        typedef DynamicAutoDiffCostFunction<BaseGridProjection<Projector>> gridProjectionCF;
        BaseGridProjection<Projector> * gridProjection;
        gridProjection = new BaseGridProjection<Projector>(projection, grid, 
                                    Transformation<double>(poseOrigBase.data()), intrinsic);
        gridProjectionCF * costFunction = new gridProjectionCF(gridProjection);
        costFunction->AddParameterBlock(6);
        costFunction->AddParameterBlock(6);
        costFunction->SetNumResiduals(2 * Nx * Ny);
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                poseCameraBase.data(), poseOrigGrid.data());
    }
    
    void estimateInitialGrid(const vector<double> & intrinsic, array<double, 6> & poseOrigGrid,
            const vector<CalibrationData> & calibDataVec)
    {
        auto const & item = calibDataVec[0];
        array<double, 6> poseCamGrid;
        estimateInitialGrid(intrinsic, item.projection, poseCamGrid);
        Transformation<double> TcamGrid(poseCamGrid.data());
        Transformation<double> TbaseCam(initialBaseCameraGuess.data());
        Transformation<double> TorigBase(item.extrinsic.data());
        Transformation<double> TorigGrid = TorigBase.compose(TbaseCam).compose(TcamGrid);
        
        TorigGrid.toArray(poseOrigGrid.data());
    }

    //TODO compare to robot calibration
    bool compute(Transformation<double> & TbaseCam, Transformation<double> & TorigGrid)
    {

        cout << "### Extrinsic parameters calibration ###" << endl;
        if (baseCalibDataVec.size() == 0)
        {
            cout << "ERROR : none of images were accepted" << endl;
            return false;
        } 
        
        array<double, 6> poseOrigGrid;
        estimateInitialGrid(intrinsicVec, poseOrigGrid, baseCalibDataVec);
        
        array<double, 6> poseCamBase = 
                Transformation<double>(initialBaseCameraGuess.data()).inverse().toArray();
        
        cout << "Global optimization" << endl;
        Problem problem;
        
        // Intrinsic init
        for (auto & item : baseCalibDataVec)
        {
            addBaseResidual(problem, intrinsicVec, item.extrinsic, item.projection,
                    poseCamBase, poseOrigGrid);
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
        
        TbaseCam = Transformation<double>(poseCamBase.data()).inverse();
        TorigGrid = Transformation<double>(poseOrigGrid.data());
    } 
    
//    bool residualAnalysis(const vector<double> & intrinsic1, 
//            const vector<double> & intrinsic2,
//            const Transformation<double> & TleftRight)
//    {
//        residualAnalysis(intrinsic1, monoCalibDataVec1);
//        residualAnalysis(intrinsic2, monoCalibDataVec2);
//        vector<CalibrationData> stereoLeftData, stereoRightData;
//        for (auto & x : stereoCalibDataVec)
//        {
//            stereoLeftData.push_back(x.getLeftMono());
//            stereoRightData.push_back(x.getRightMono(TleftRight));
//        }
//        residualAnalysis(intrinsic1, stereoLeftData);
//        residualAnalysis(intrinsic2, stereoRightData);
//    }
};

