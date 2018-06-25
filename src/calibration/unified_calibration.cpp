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

#include "calibration/unified_calibration.h"

#include "std.h"
#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "ceres.h"
#include "json.h"
#include "except.h"
#include "timer.h"

#include <glog/logging.h>

#include "calibration/calib_cost_functions.h"
#include "calibration/odometry_cost_function.h"
#include "calibration/corner_detector.h"
#include "projection/generic_camera.h"
#include "projection/eucm.h"
#include "projection/ucm.h"
#include "projection/mei.h"

bool GenericCameraCalibration::compute()
{
    //run the solver
    Solver::Options options;
//    options.check_gradients = true;
//    options.gradient_check_relative_precision = 0.01;
//    options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
    options.max_num_iterations = 1000;
    options.function_tolerance = 1e-15;
    options.gradient_tolerance = 1e-15;
    options.parameter_tolerance = 1e-15;
//    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &globalProblem, &summary);
    cout << summary.FullReport() << endl;
    
    cout << "Intrinsic parameters :" << endl;
    for (auto & x : intrinsicMap)
    {
        cout << x.first << " : ";
        for (int i = 0; i < cameraMap[x.first]->numParams(); i++)
        {
            cout << x.second[i] << "  ";
        }
        cout << endl;
        
    }
    
    cout << "Local extrinsic parameters :" << endl;
    for (auto & sequenceMap : sequenceTransformMap)
    {
        cout << "Sequence : " << sequenceMap.first << endl;
        int i = 0;
        for (auto & x : sequenceMap.second)
        {
            cout << i++ << " : " << Transf(x.data()) << endl;
        }
    }
    
    cout << "Global extrinsic parameters :" << endl;
    for (auto & x : globalTransformMap)
    {
        cout << x.first << " : " << Transf(x.second.data()) << endl;
    }
    
    for (int i = 0; i < dataVec.size(); i++)
    {
        writeImageResidual(dataVec[i], "image_error_" + to_string(i) + ".txt");
    }
}

void GenericCameraCalibration::parseTransforms()
{
    for (auto & transInfo : root.get_child("transformations"))
    {
        const string name = transInfo.second.get<string>("name");
        transformInfoMap[name] = TransformInfo();
        auto & info = transformInfoMap[name];
        info.global = transInfo.second.get<bool>("global");
        info.prior = transInfo.second.get<bool>("prior");
        info.constant = transInfo.second.get<bool>("constant");
        if (info.constant and not info.prior)
        {
            throw runtime_error(name + " is constant but there is no prior");
        }
        
        info.initialized = false;
        
        if (info.global)
        {
            globalTransformMap[name] = Array6d();
            if (info.prior)
            {
                auto transf = readTransform(transInfo.second.get_child("value"));
                transf.toArray(globalTransformMap[name].data());
            }
        }
        else 
        {
            sequenceTransformMap[name] = vector<Array6d>();
            sequenceInitMap[name] = vector<bool>();
            if (info.prior)
            {
                auto & valVec = sequenceTransformMap[name];
                for (auto & val : transInfo.second.get_child("value"))
                {
                    auto transf = readTransform(val.second);
                    valVec.emplace_back(transf.toArray());
                }
            }
        }
    }
}

void GenericCameraCalibration::parseCameras()
{
    for (auto & cameraInfo : root.get_child("cameras"))
    {
        const string name = cameraInfo.second.get<string>("name");
        cameraConstantMap[name] = cameraInfo.second.get<bool>("constant");
        intrinsicMap[name] = vector<double>();
        auto & intrinsicVec = intrinsicMap[name];
        const string cameraType = cameraInfo.second.get<string>("type");
        for (auto & x : cameraInfo.second.get_child("value"))
        {
            intrinsicVec.push_back(x.second.get_value<double>());
        }
        //TODO add other projection models in the future
        if (cameraType == "eucm")
        {
            cout << "Model : EUCM" << endl;
            if (intrinsicVec.size() != 6)
            {
                throw runtime_error("invalid number of intrinsic parameters");
            }
            cameraMap[name] = new EnhancedCamera(intrinsicVec.data());
        }
        else if (cameraType == "ucm")
        {
            cout << "Model : UCM" << endl;
            if (intrinsicVec.size() != 5)
            {
                throw runtime_error("invalid number of intrinsic parameters");
            }
            cameraMap[name] = new UnifiedCamera(intrinsicVec.data());
        }
        else if (cameraType == "mei")
        {
            cout << "Model : MEI" << endl;
            if (intrinsicVec.size() != 10)
            {
                throw runtime_error("invalid number of intrinsic parameters");
            }
            cameraMap[name] = new MeiCamera(intrinsicVec.data());
        }
        else
        {
            throw runtime_error("invalid camera model name");
        }
    }
}

void GenericCameraCalibration::initTransformChainInfo(ImageData & data, const ptree & node)
{
    data.cameraName = node.get<string>("camera");
    for (auto & flag : node.get_child("parameters"))
    {
        string flagName = flag.second.get_value<string>();
        if (flagName == "check_extraction") data.checkExtraction = true;
        else if (flagName == "improve_detection") data.improveDetection = true;
        else if (flagName == "show_outliers") data.showOutliers = true;
        else if (flagName == "user_guided") data.userGuided = true;
        else if (flagName == "do_not_solve") data.doNotSolve = true;
        else if (flagName == "do_not_solve_global") data.doNotSolveGlobal = true;
        else if (flagName == "save_outlire_images") data.saveOutlierImages = true;
        else if (flagName == "draw_improved") 
        {
            data.drawImproved = true;
            data.drawScale = node.get<double>("draw_scale");
        }
        else
        {
            cout << "WARNING : UNKNOWN FLAG -- " << flagName << endl;
        }
    }
    cout <<"Camera : " <<  data.cameraName << endl;
    cout <<"Transformations : ";
    for (auto & transInfo : node.get_child("transform_chain"))
    {
        data.transNameVec.push_back(transInfo.second.get<string>("name"));
        cout << data.transNameVec.back();
        if (transInfo.second.get<bool>("direct")) 
        {
            data.transStatusVec.push_back(TRANSFORM_DIRECT);
        }
        else
        {
            data.transStatusVec.push_back(TRANSFORM_INVERSE);
            cout << "_inv";
        }
        cout <<  "   ";
    }
    //count transformation sequences
    int sequenceCount = 0;
    for (auto & name : data.transNameVec)
    {
        if (not transformInfoMap[name].global) sequenceCount++;
    }
    if (sequenceCount != 1) throw runtime_error("not one sequences in a transform chain");
    
    cout << endl;
}


void GenericCameraCalibration::initGridIR(ImageData & data, const ptree & node)
{   
    data.Nx = 2;
    data.Ny = 2;
    data.sqSize = -1;
    data.board.clear();
    data.board.reserve(data.Nx * data.Ny);
    data.idxUL = node.get<int>("object.corner_ul");
    data.idxUR = node.get<int>("object.corner_ur");
    data.idxBL = node.get<int>("object.corner_bl");
    data.idxBR = node.get<int>("object.corner_br");
    for (auto & x : node.get_child("object.points"))
    {
        vector<double> pt = readVector<double>(x.second);
        data.board.emplace_back(pt[0], pt[1], pt[2]);
    }
}

void GenericCameraCalibration::readCorners(ImageData & data, const ptree & node)
{
    data.useImages = false;
    data.imageWidth = node.get<int>("image_width");
    data.imageHeight = node.get<int>("image_height");
    ptree dataFile;
    read_json(node.get<string>("data_file"), dataFile);
    string cameraID = node.get<string>("camera");
    for (auto & dataPoint : dataFile)
    {
        data.detectedCornersVec.emplace_back();
        Vector2dVec & cornerVec = data.detectedCornersVec.back();
        for (auto & x : dataPoint.second)
        {
            if (x.second.get<string>("camera") == cameraID)
            {
                for (auto & y : x.second.get_child("points"))
                {
                    vector<double> pt = readVector<double>(y.second);
                    cornerVec.emplace_back(pt[0], pt[1]);
                }
                break;
            }
        }
    }
}

void GenericCameraCalibration::initGrid(ImageData & data, const ptree & node)
{   
    data.Nx = node.get<int>("object.cols");
    data.Ny = node.get<int>("object.rows");
    data.sqSize = node.get<double>("object.size");
    data.board.clear();
    data.board.reserve(data.Nx * data.Ny);
    for (int i = 0; i < data.Ny; i++)
    {
        for (int j = 0; j < data.Nx; j++)
        {
           data.board.emplace_back(data.sqSize * j, data.sqSize * i, 0); 
        }
    }
    data.idxUL = 0;
    data.idxUR = data.Nx - 1;
    data.idxBL = data.Nx * (data.Ny - 1);
    data.idxBR = data.Nx * data.Ny - 1;
    
    //fill up detectedCornersVec which stores all the extracted grids
    data.useImages = true;
    const string prefix = node.get<string>("images.prefix");
    data.detectedCornersVec.clear();
    for (auto & x : node.get_child("images.names"))
    {
        const string filename = x.second.get_value<string>();
        data.imageNameVec.emplace_back(prefix + filename);
    }
    extractGridProjections(data);

}

Transf GenericCameraCalibration::getInitTransform(Transf xi,
            const string & initName, const ImageData & data, const int transfIdx)
{
    for (int i = 0; i < data.transNameVec.size(); i++)
    {
        const string & name = data.transNameVec[i];
        if (name == initName) break;
        else if (data.transStatusVec[i] == TRANSFORM_DIRECT)
        {
            xi = getTransform(name, transfIdx).inverseCompose(xi);
        }
        else if (data.transStatusVec[i] == TRANSFORM_INVERSE)
        {
            xi = getTransform(name, transfIdx).compose(xi);
        }
    }
    for (int i = data.transNameVec.size() - 1; i >= 0; i--)
    {
        const string & name = data.transNameVec[i];
        if (name == initName)
        {
            if (data.transStatusVec[i] == TRANSFORM_INVERSE)
            {
                xi = xi.inverse();
            }
            break;
        }
        else if (data.transStatusVec[i] == TRANSFORM_DIRECT)
        {
            xi = xi.composeInverse(getTransform(name, transfIdx));
        }
        else if (data.transStatusVec[i] == TRANSFORM_INVERSE)
        {
            xi = xi.compose(getTransform(name, transfIdx));
        }
    }
    return xi;
}

bool GenericCameraCalibration::addResiduals(const string & infoFileName)
{
    read_json(infoFileName, root);
    parseTransforms();
    parseCameras();
    parseData();
}

void GenericCameraCalibration::initGlobalTransform(const ImageData & data, const string & name)
{
    Problem problem;
    for (int transfIdx = 0; transfIdx < data.detectedCornersVec.size(); transfIdx++)
    {
        if (data.detectedCornersVec[transfIdx].empty()) continue;
        
        vector<double*> ptrVec;
        for (int i = 0; i < data.transNameVec.size(); i++)
        {
            const string & name = data.transNameVec[i];
            ptrVec.push_back(getTransformData(name, transfIdx).data());
        }

        //add a residual
        GenericProjectionJac * costFunction = new GenericProjectionJac(data.detectedCornersVec[transfIdx],
                    data.board, cameraMap[data.cameraName], data.transStatusVec);

        switch (ptrVec.size())
        {
        case 0:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data());
            break;
        case 1:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data(), ptrVec[0]); //FIXME
            break;
        case 2:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data(), ptrVec[0], ptrVec[1]);
            break;
        case 3:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data(), ptrVec[0], ptrVec[1], ptrVec[2]);
            break;
        case 4:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data(),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3]);
            break;
        case 5:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    intrinsicMap[data.cameraName].data(),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], ptrVec[4]);
            break;
        }
        
        //set everything to constant except for the transform to be initialize
        for (int i = 0; i < data.transNameVec.size(); i++)
        {
            if (data.transNameVec[i] != name)
            {
                problem.SetParameterBlockConstant(ptrVec[i]);
            }
        }
    }
    //initrinsics are constant as well
    problem.SetParameterBlockConstant(intrinsicMap[data.cameraName].data());
    
    //solve
    Solver::Options options;
//    options.check_gradients = true;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 500;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
//    cout << summary.FullReport() << endl;
//    cout << getTransform(name) << endl;
}

void GenericCameraCalibration::initTransforms(const ImageData & data, const string & initName)
{
    if (initName == "none") return;
    else if (transformInfoMap.find(initName) == transformInfoMap.end())
    {
        throw runtime_error(initName + " does not exist, impossible to initialize");
    }
    
    auto nameIter = find(data.transNameVec.begin(), data.transNameVec.end(), initName);
    if (nameIter == data.transNameVec.end())
    {
        throw runtime_error(initName + " does not belong to the transform chain");
    }
    
    if  (transformInfoMap[initName].prior) 
    {
        throw runtime_error(initName + " has a prior value");
    }
    
//    if  (transformInfoMap[initName].initialized) 
//    {
//        throw runtime_error(initName + " has already been initialized");
//    }
    
    
//    if (not transformInfoMap[initName].global and sequenceTransformMap[initName].size() != 0)
//    {
//        throw runtime_error(initName + " has already been initialized");
//    }
    
    transformInfoMap[initName].initialized = true;
    //FIXME if the transformation is partially initialized then it breaks !!!!! BUG
    for (auto & x : data.transNameVec)
    {
        if ( not (transformInfoMap[x].prior xor transformInfoMap[x].initialized) )
        {
            const string errMsg = " is not initialized. Cannot initialize more than one transform at a time";
            throw runtime_error(x + errMsg);
        }
    }
    
    //do the initialization
    if (not transformInfoMap[initName].global)
    {
        const bool IS_ALLOCATED = (sequenceTransformMap[initName].size() != 0);
        if (not IS_ALLOCATED)
        {
            sequenceTransformMap[initName].reserve(data.detectedCornersVec.size());
        }
        
        for (int transfIdx = 0; transfIdx < data.detectedCornersVec.size(); transfIdx++)
        {
            if (not IS_ALLOCATED)
            {
                sequenceTransformMap[initName].push_back(Array6d{0, 0, 1, 0, 0, 0});  
                sequenceInitMap[initName].push_back(false);
            }
            
            if (not data.detectedCornersVec[transfIdx].empty() and not sequenceInitMap[initName][transfIdx])
            {
//                cout << "WARNING : " << initName << " " << transfIdx
//                     << " is not initialized, no board extracted" << endl;
                auto xi = estimateInitialGrid(data, transfIdx);
                xi = getInitTransform(xi, initName, data, transfIdx);
                sequenceTransformMap[initName][transfIdx] = xi.toArray();
                sequenceInitMap[initName][transfIdx] = true;
            }
        }
    }
    else
    {
        int transfIdx = data.getFirstExtractedIdx();
        auto xi = estimateInitialGrid(data, transfIdx);
        cout << "INITI VALUE IN CAMERA FRAME " << endl;
        cout << xi << endl;
        xi = getInitTransform(xi, initName, data, transfIdx);
        cout << "INITI TRANSFORM " << endl;
        cout << xi << endl;
        xi.toArray(globalTransformMap[initName].data());
        if (data.detectedCornersVec.size() > 1 and not data.doNotSolve) initGlobalTransform(data, initName);
    }
}

void GenericCameraCalibration::addGridResidualBlocks(const ImageData & data)
{
    if (data.doNotSolveGlobal) return;

    for (int transfIdx = 0; transfIdx < data.detectedCornersVec.size(); transfIdx++)
    {
        if (data.detectedCornersVec[transfIdx].empty()) continue;

        // make the vector of pointers to the transformation data
        vector<double*> ptrVec;
        for (int i = 0; i < data.transNameVec.size(); i++)
        {
            const string & name = data.transNameVec[i];
            ptrVec.push_back(getTransformData(name, transfIdx).data());
        }

        //TODO create header calibration_data_structures.h
        //add a residual
        GenericProjectionJac * costFunction = new GenericProjectionJac(data.detectedCornersVec[transfIdx],
                    data.board, cameraMap[data.cameraName], data.transStatusVec);
        double * intrinsicPtr = intrinsicMap[data.cameraName].data();
        
        switch (ptrVec.size())
        {
        case 0:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr);
            break;
        case 1:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr, ptrVec[0]); //FIXME
            break;
        case 2:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr, ptrVec[0], ptrVec[1]);
            break;
        case 3:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr, ptrVec[0], ptrVec[1], ptrVec[2]);
            break;
        case 4:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr,
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3]);
            break;
        case 5:
            globalProblem.AddResidualBlock(costFunction, NULL, //new SoftLOneLoss(1),
                    intrinsicPtr,
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], ptrVec[4]);
            break;
        default:
            throw runtime_error("the transform chain is too long (5 transforms at max are supproted)");
        }
        
//        
//        switch (ptrVec.size())
//        {
//        case 0:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2), intrinsicPtr);
//            break;
//        case 1:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2),
////                    ptrVec[0], intrinsicPtr); //FIXME
//                    intrinsicPtr, ptrVec[0]);
//            break;
//        case 2:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2),
//                    ptrVec[0], ptrVec[1], intrinsicPtr);
//            break;
//        case 3:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2),
//                    ptrVec[0], ptrVec[1], ptrVec[2], intrinsicPtr);
//            break;
//        case 4:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2),
//                    ptrVec[0], ptrVec[1], ptrVec[2],
//                    ptrVec[3], intrinsicPtr);
//            break;
//        case 5:
//            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(2),
//                    ptrVec[0], ptrVec[1], ptrVec[2],
//                    ptrVec[3], ptrVec[4], intrinsicPtr);
//            break;
//        default:
//            throw runtime_error("the transform chain is too long (5 transforms at max are supproted)");
//        }
        
        //set constant transformations
        for (int i = 0; i < data.transNameVec.size(); i++)
        {
            const string & name = data.transNameVec[i];
            if (transformInfoMap[name].constant)
            {
                globalProblem.SetParameterBlockConstant(ptrVec[i]);
            }
            ptrVec.push_back(getTransformData(name, transfIdx).data());
        }

        if (cameraConstantMap[data.cameraName]) 
        {
            globalProblem.SetParameterBlockConstant(intrinsicMap[data.cameraName].data());
        }
        else
        {        
            //set the limits on intrinsics
            for (int i = 0; i < intrinsicMap[data.cameraName].size(); i++)
            {
                globalProblem.SetParameterLowerBound(intrinsicMap[data.cameraName].data(), i,
                            cameraMap[data.cameraName]->lowerBound(i));
                globalProblem.SetParameterUpperBound(intrinsicMap[data.cameraName].data(), i,
                            cameraMap[data.cameraName]->upperBound(i));
            }
        }
    }
}
    
void GenericCameraCalibration::parseData()
{
    for (auto & dataInfo : root.get_child("data"))
    {
        const string dataType = dataInfo.second.get<string>("type");
        if (dataType == "images")
        {
            //load calibration data
            dataVec.emplace_back();
            initTransformChainInfo(dataVec.back(), dataInfo.second);
            initGrid(dataVec.back(), dataInfo.second);
            
            //init variables and add residuals to the problem
            initTransforms(dataVec.back(), dataInfo.second.get<string>("init"));
            addGridResidualBlocks(dataVec.back());
        }
        if (dataType == "ir_data")
        {
            //load calibration data
            dataVec.emplace_back();
            initTransformChainInfo(dataVec.back(), dataInfo.second);
            initGridIR(dataVec.back(), dataInfo.second);
            
            
            readCorners(dataVec.back(), dataInfo.second);
            //init variables and add residuals to the problem
            initTransforms(dataVec.back(), dataInfo.second.get<string>("init"));
            addGridResidualBlocks(dataVec.back());
        }
        if (dataType == "odometry_intrinsic")
        {
            const string transformName = dataInfo.second.get<string>("transform");
            if (transformInfoMap.find(transformName) == transformInfoMap.end())
            {
                throw runtime_error(transformName + " has not been declared");
            }
            if (transformInfoMap[transformName].global)
            {
                throw runtime_error(transformName + " is global. Odometry must be a sequence");
            }
            
            const double errV = dataInfo.second.get<double>("err_v"); //relative error in speed
            const double errW = dataInfo.second.get<double>("err_w"); //relative error in rotation
            const double lambda = dataInfo.second.get<double>("lambda"); //relative error in rotation
            
            intrinsicMap[transformName] = vector<double>();
            intrinsicMap[transformName].push_back(dataInfo.second.get<double>("radius_left"));
            intrinsicMap[transformName].push_back(dataInfo.second.get<double>("radius_right"));
            intrinsicMap[transformName].push_back(dataInfo.second.get<double>("track_gauge"));
            vector<Vector2dVec> deltaQVec;
            
            //read the odometry increment measurements
            ptree dataFile;
            read_json(dataInfo.second.get<string>("data_file"), dataFile);
            for (auto & dataPoint : dataFile)
            {
                deltaQVec.emplace_back();
                for (auto & x : dataPoint.second)
                {
                    vector<double> pt = readVector<double>(x.second);
                    deltaQVec.back().emplace_back(pt[0], pt[1]);
                }
            }
            
            //use the odometry as initial values
            if (dataInfo.second.get<bool>("init"))
            {
                cout << transformName << endl;
//                for (auto & xx : sequenceTransformMap[transformName])
//                {
//                    for (auto & x : xx)
//                    {
//                        cout << x << "   ";
//                    }
//                    cout << endl;
//                }
                if (not sequenceTransformMap[transformName].empty())
                {
                    throw runtime_error(transformName + " has already been initialized");
                }
                transformInfoMap[transformName].initialized = true;
                sequenceTransformMap[transformName].emplace_back(array<double, 6>{0, 0, 0, 0, 0, 0});
                sequenceInitMap[transformName].push_back(true);
                //The init will be done at the next step
            }
            
            //add the cost functions
            for (int i = 0; i < deltaQVec.size(); i++)
            {
                OdometryCost * costFunction = new OdometryCost(errV, errW, lambda,
                                                    deltaQVec[i], intrinsicMap[transformName].data());
                
                if (dataInfo.second.get<bool>("init"))
                {
                    Transf xiPrev(sequenceTransformMap[transformName].back().data());
                    Transf xiPrior = xiPrev.compose(costFunction->_zetaPrior);
                    sequenceTransformMap[transformName].emplace_back(xiPrior.toArray());
                    sequenceInitMap[transformName].push_back(true);
                }   
                    
                globalProblem.AddResidualBlock(costFunction, NULL,
                    sequenceTransformMap[transformName][i].data(),
                    sequenceTransformMap[transformName][i + 1].data(),
                    intrinsicMap[transformName].data()); 
                
            }
            if (dataInfo.second.get<bool>("anchor"))
            {
                globalProblem.SetParameterBlockConstant(sequenceTransformMap[transformName][0].data());
            }
        }
        else if (dataType == "odometry")
        {
            const string transformName = dataInfo.second.get<string>("transform");
            if (transformInfoMap.find(transformName) == transformInfoMap.end())
            {
                throw runtime_error(transformName + " has not been declared");
            }
            if (transformInfoMap[transformName].global)
            {
                throw runtime_error(transformName + " is global. Odometry must be a sequence");
            }
            
            const double errV = dataInfo.second.get<double>("err_v"); //relative error in speed
            const double errW = dataInfo.second.get<double>("err_w"); //relative error in rotation
            const double lambda = dataInfo.second.get<double>("lambda"); //relative error in rotation
            //read out the transformations
            vector<Transf> odometryVec;
            for (auto & odomItem : dataInfo.second.get_child("value"))
            {
                odometryVec.emplace_back(readTransform(odomItem.second));
                auto & xi = odometryVec.back();
//                xi.rot() += 0.002 * Vector3d::Random();
//                xi.trans() += 0.01 * Vector3d::Random();
            }
            
            //use the odometry as initial values
            if (dataInfo.second.get<bool>("init"))
            {
                cout << transformName << endl;
//                for (auto & xx : sequenceTransformMap[transformName])
//                {
//                    for (auto & x : xx)
//                    {
//                        cout << x << "   ";
//                    }
//                    cout << endl;
//                }
                if (not sequenceTransformMap[transformName].empty())
                {
                    throw runtime_error(transformName + " has already been initialized");
                }
                transformInfoMap[transformName].initialized = true;
                for (auto & xi : odometryVec)
                {
                    sequenceTransformMap[transformName].emplace_back(xi.toArray());
                    sequenceInitMap[transformName].push_back(true);
                }
            }
            
            //add the cost functions
            for (int i = 0; i < odometryVec.size() - 1; i++)
            {
                OdometryPrior * costFunction = new OdometryPrior(errV, errW, lambda,
                                odometryVec[i], odometryVec[i + 1]);
                                
                globalProblem.AddResidualBlock(costFunction, NULL,
                    sequenceTransformMap[transformName][i].data(),
                    sequenceTransformMap[transformName][i + 1].data()); 
                
            }
            if (dataInfo.second.get<bool>("anchor"))
            {
                globalProblem.SetParameterBlockConstant(sequenceTransformMap[transformName][0].data());
            }
        }
        else if (dataType == "transformation_prior")
        {
            const string transformName = dataInfo.second.get<string>("transform");
            
            if (transformInfoMap.find(transformName) == transformInfoMap.end())
            {
                throw runtime_error(transformName + " has not been declared");
            }
            if (not transformInfoMap[transformName].prior)
            {
                throw runtime_error(transformName + " must have a prior value");
            }
            vector<double> stiffnessVec;
            for (auto & x : dataInfo.second.get_child("stiffness"))
            {
                stiffnessVec.push_back(x.second.get_value<double>());
            }
            
            double * transformData = getTransformData(transformName).data();
            TransformationPrior * costFunction = new TransformationPrior(stiffnessVec.data(), transformData);
            globalProblem.AddResidualBlock(costFunction, NULL, transformData); 
        }
    }
}

//FIXME
vector<int> xVec;
vector<int> yVec;
void onMouse( int event, int x, int y, int, void* )
{
    if( event != CV_EVENT_LBUTTONDOWN )
        return;

    xVec.push_back(x);
    yVec.push_back(y);
}

// OLD VECSION FOR OPENCV COMPARISON --- segfault 0_o
/*
void GenericCameraCalibration::extractGridProjections(ImageData & data)
{
    Timer timer;
    const int flags = cv::CALIB_CB_ADAPTIVE_THRESH;
    Size patternSize(data.Nx, data.Ny);
    string sequenceName;
    for (auto & name : data.transNameVec)
    {
        if (transformInfoMap[name].global == false)
        {
            sequenceName = name;
            break;
        }
    }
    bool initialized = transformInfoMap[sequenceName].initialized;
    const vector<bool> & initVec = sequenceInitMap[sequenceName];
    int countSuccess = 0;
    for (int i = 0; i < data.imageNameVec.size(); i++)
    {
        cout << "." << flush;
        data.detectedCornersVec.emplace_back();
        
        //grid has not been found on the corresponding image in the ini sequence
        if (initialized and not initVec[i]) continue;
        
        const string & fileName = data.imageNameVec[i];
        Mat8u frame = imread(fileName, 0);
        if (frame.empty())
        {
            cout << fileName << " : ERROR, file not found" << endl;
            continue;
        }
        vector<cv::Point2f> centers;
        bool patternIsFound = findChessboardCorners(frame, patternSize, centers, flags);
        if (not patternIsFound)
        {
            //USER-guided detection
            if (data.userGuided)
            {
                xVec.clear();
                yVec.clear();
                
                cout << "CB detection failed" << endl;
                cout << "select UL, then BR corners of the area with CB" << endl;

                imshow("Select", frame);
                setMouseCallback("Select", onMouse);
                while ( xVec.size() < 2) waitKey(50);

                Mat8u subframe = frame.colRange(xVec[0], xVec[1]).rowRange(yVec[0], yVec[1]);
                double ratio = 1;
                const double ROWS_MIN = 120;
                const double COLS_MIN = 160;
                while (subframe.rows * ratio < ROWS_MIN or subframe.cols * ratio < COLS_MIN)
                {
                    ratio += 1;
                }
                if (ratio != 1) resize(subframe, subframe, Size(0, 0), ratio, ratio);
                
                patternIsFound = findChessboardCorners(subframe, patternSize,
                                        centers, flags);
                imshow("Select", subframe);
                waitKey();
                if (not patternIsFound)
                {
                    cout << fileName << " : ERROR, pattern not found" << endl;
                    continue;
                }
                for (auto & pt : centers)
                {
                    pt.x = pt.x / ratio + xVec[0];
                    pt.y = pt.y / ratio + yVec[0];
                }
            }
            else
            {
                cout << fileName << " : ERROR, pattern not found" << endl;
                continue;
            }
        }
        
        if (data.checkExtraction)
        {
            //FIXME 
//            F0310 19:01:41.823675  1227 cubic_interpolation.h:72] Check failed: x >= 0.0 (-nan vs. 0) 
// *** Check failure stack trace: ***
// Aborted (core dumped)
            Mat8uc3 cornerImg;
            cvtColor(frame, cornerImg, CV_GRAY2BGR);
            drawChessboardCorners(cornerImg, patternSize, Mat(centers), patternIsFound);
            imshow("corners", cornerImg);
            char key = waitKey();
            if (key == 'n' or key == 'N')
            {
                cout << fileName << " : ERROR, pattern not accepted" << endl;
                continue;
            }
        }
        
        auto & cornerVec = data.detectedCornersVec.back();
        cornerVec.reserve(centers.size());
        
        Mat8uc3 cornerImg;
        const double & K = data.drawScale;
        if (data.drawImproved)
        {
            cvtColor(frame, cornerImg, CV_GRAY2BGR);
            resize(cornerImg, cornerImg, Size(0, 0), K, K);
        }
        for (auto pt : centers)
        {
            cornerVec.emplace_back(pt.x, pt.y);
            if (data.drawImproved) 
            {
                cross(cornerImg, pt.x * K + 0.5*K, pt.y * K + 0.5*K, 25, 255, 3);
            }
        }
        countSuccess++;
        
//        if (data.improveDetection)
//        {
//            double minDist = findMinDistance(cornerVec, data.Ny, data.Nx);
//            CornerDetector detector(data.Nx, data.Ny, 3, true, false);
//            detector.improveCorners(cornerVec);
//        }
        
        if (data.drawImproved) 
        {
            for (auto pt : cornerVec)
            {
                cross(cornerImg, pt[0] * K + 0.5*K, pt[1] * K + 0.5*K, 25, Scalar(0, 0, 255), 3);
            }
            imshow("corners", cornerImg);
            waitKey();
        }
    }
    double telapsed = timer.elapsed();
    cout << endl;
    cout << "DETECTION RATE : " << countSuccess << " of " 
            << data.imageNameVec.size() << " detected" << endl;
    cout << "ELAPSED : " << telapsed << "      or per image : " << telapsed / data.imageNameVec.size() << endl;
}*/



void GenericCameraCalibration::extractGridProjections(ImageData & data)
{
    Timer timer;
    CornerDetector detector(data.Nx, data.Ny, 3, data.improveDetection);
    
    string sequenceName;
    for (auto & name : data.transNameVec)
    {
        if (transformInfoMap[name].global == false)
        {
            sequenceName = name;
            break;
        }
    }
    bool initialized = transformInfoMap[sequenceName].initialized;
    const vector<bool> & initVec = sequenceInitMap[sequenceName];
    int countSuccess = 0;
    for (int i = 0; i < data.imageNameVec.size(); i++)
    {
        cout <<data.imageNameVec[i] << endl;
//        cout << "." << flush;
        data.detectedCornersVec.emplace_back();
        
        //grid has not been found on the corresponding image in the ini sequence
        
        const string & fileName = data.imageNameVec[i];
        
        if (initialized and not initVec[i])
        {
            cout << fileName << " : ERROR, the pattern has not been found on the corresponding image" << endl;
            continue;
        }
        
        Mat8u frame = imread(fileName, 0);
        if (frame.empty())
        {
            cout << fileName << " : ERROR, file not found" << endl;
            continue;
        }
        detector.setImage(frame);
        Vector2dVec patternVec;
        bool patternIsFound = detector.detectPattern(patternVec);
        if (not patternIsFound)
        {
            cout << fileName << " : ERROR, pattern not found" << endl;
            continue;
        }
        
        if (data.checkExtraction)
        {
            Mat8u cornerImg;
            frame.copyTo(cornerImg);
            
            drawPoints(cornerImg, patternVec);
            imshow("corners", cornerImg);
            char key = waitKey();
            if (key == 'n' or key == 'N')
            {
                cout << fileName << " : ERROR, pattern not accepted" << endl;
                continue;
            }
        }
        
        countSuccess++;
        data.detectedCornersVec.back() = patternVec;
    }
    double telapsed = timer.elapsed();
    cout << endl;
    cout << "DETECTION RATE : " << countSuccess << " of " 
            << data.imageNameVec.size() << " detected" << endl;
    cout << "ELAPSED : " << telapsed << "      or per image : " << telapsed / data.imageNameVec.size() << endl;
}


Transf GenericCameraCalibration::estimateInitialGrid(const ImageData & data, const int gridIdx)
{
    auto & cornerVec = data.detectedCornersVec[gridIdx];
    
    ICamera * cam = cameraMap[data.cameraName];
    
    
    array<double, 6> xiArr{0, 0, 1, 0, 0, 0};
    
    
    Vector2d ptUL = cornerVec[data.idxUL];
    Vector2d ptUR = cornerVec[data.idxUR];
    Vector2d ptBL = cornerVec[data.idxBL];
    Vector2d ptBR  = cornerVec[data.idxBR];
    Vector3d XUL, XUR, XBL, XBR;
    
    cam->reconstructPoint(ptUL, XUL);
    cam->reconstructPoint(ptUR, XUR);
    cam->reconstructPoint(ptBL, XBL);
    cam->reconstructPoint(ptBR, XBR);
    
    XUL.normalize();
    XUR.normalize();
    XBR.normalize();
    XBL.normalize();
    
    //initial position
    Vector3d exU(XUR - XUL);
    Vector3d exB(XBR - XBL);
    Vector3d eyL(XBL - XUL);
    Vector3d eyR(XBR - XUR);
    double exModelU = (data.board[data.idxUR] - data.board[data.idxUL]).norm();
    double exModelB = (data.board[data.idxBR] - data.board[data.idxBL]).norm();
    double eyModelL = (data.board[data.idxBL] - data.board[data.idxUL]).norm();
    double eyModelR = (data.board[data.idxBR] - data.board[data.idxUR]).norm(); 
    
    double scaleXU = exModelU / exU.norm();
    double scaleXB = exModelB / exB.norm();
    double scaleYL = eyModelL / eyL.norm();
    double scaleYR = eyModelR / eyR.norm();
    
    Vector3d pos = XUL * min(scaleXU, scaleYL);
    copy(pos.data(), pos.data() + 3, xiArr.data());
    
    //initial orientation
    //compute the board basis and make a rotation matrix
    /*
    Vector3d ez = (X1 + X2) * 0.5; // ex.dot(ez) = 0
    ex1.normalize();
    ez.normalize();
    Vector3d ey = ez.cross(ex);
    Matrix3d R;
    R << ex1, ey, ez;
    */
    Vector3d posx = XUR * min(scaleXU, scaleYR);
    Vector3d posy = XBL * min(scaleXB, scaleYL);
    Vector3d ex = posx - pos;
    Vector3d ey = posy - pos;
    ex.normalize();
    //make ey perpendicular
    ey = (Matrix3d::Identity() - ex * ex.transpose()) * ey;
    ey.normalize();
    
    Vector3d ez = ex.cross(ey);
    Matrix3d R;
    R << ex, ey, ez;
    
    Vector3d rot = rotationVector(R);
    
    copy(rot.data(), rot.data() + 3, xiArr.data() + 3);
    
    if (not data.doNotSolve)
    {
        Problem problem;
        GenericProjectionJac * costFunction = new GenericProjectionJac(cornerVec, data.board,
                cam, {TRANSFORM_DIRECT});
                
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(25),
                intrinsicMap[data.cameraName].data(), xiArr.data());
        problem.SetParameterBlockConstant(intrinsicMap[data.cameraName].data());
        Solver::Options options;
//        options.check_gradients = true;
//        options.gradient_check_relative_precision = 1;
        options.max_num_iterations = 500;
        Solver::Summary summary;
        
        Solve(options, &problem, &summary);
//        cout << "IMAGE NAME : " << data.imageNameVec[gridIdx] << endl;
//        cout << summary.FullReport() << endl;
    }
    
    return Transf(xiArr.data());
}

void GenericCameraCalibration::computeTransforms(const ImageData & data, vector<Transf> & transfVec) const
{
    transfVec.clear();
    transfVec.reserve(data.detectedCornersVec.size());
    for (int transfIdx = 0; transfIdx < data.detectedCornersVec.size(); transfIdx++)
    {
        transfVec.emplace_back(0, 0, 0, 0, 0, 0);
        auto & xi = transfVec.back();
        for (int i = 0; i < data.transNameVec.size(); i++)
        {
            const string & name = data.transNameVec[i];
            if (data.transStatusVec[i] == TRANSFORM_DIRECT)
            {
                xi = xi.compose(getTransform(name, transfIdx));
//                cout << getTransform(name, transfIdx) << endl;
            }
            else if (data.transStatusVec[i] == TRANSFORM_INVERSE)
            {
                xi = xi.composeInverse(getTransform(name, transfIdx));
//                cout << "inv " <<  getTransform(name, transfIdx) << endl;
            }
        }
    }
}


void GenericCameraCalibration::writeImageResidual(const ImageData & data, const string & fileName) const
{
    vector<Transf> transfVec;
    computeTransforms(data, transfVec);
    ofstream residualFile(fileName);
    
    //delete cam to be called later
    ICamera * cam = cameraMap.find(data.cameraName)->second->clone();
    cam->setParameters(intrinsicMap.find(data.cameraName)->second.data());
    
    for (int transfIdx = 0; transfIdx < transfVec.size(); transfIdx++)
    {
        if (data.detectedCornersVec[transfIdx].empty()) continue;
        
        Vector3dVec boardCam;
        transfVec[transfIdx].transform(data.board, boardCam);
        Vector2dVec projectedVec;
        //TODO projectPointCloud to project
        cam->projectPointCloud(boardCam, projectedVec);
        
        double stdAcc = 0;
        Vector2dVec inlierVec, outlierVec;
        vector<double> errVec;
        for (int i = 0; i < projectedVec.size(); i++)
        {
            Vector2d err = data.detectedCornersVec[transfIdx][i] - projectedVec[i];
            residualFile << err.transpose() << "   " << projectedVec[i].transpose() 
                        << "   " << transfVec[transfIdx] << endl;
            stdAcc += err.squaredNorm();
        }
        
        double sigma = sqrt(stdAcc / (projectedVec.size() - 2));
        for (int i = 0; i < projectedVec.size(); i++)
        {
            Vector2d err = data.detectedCornersVec[transfIdx][i] - projectedVec[i];
            double errNorm = err.norm();
//            if (errNorm < 3.6 * sigma and errNorm < 1.) // px is for the case when all the points are off
//            {                                        // we are looking for the sub-pixel precision
//                inlierVec.push_back(projectedVec[i]);
//            }
//            else
            {
                outlierVec.push_back(projectedVec[i]);
                errVec.push_back(err.norm());
            }
        }
        
        if ((data.showOutliers or data.saveOutlierImages) and not outlierVec.empty())
        {
            cout << "Sample #" << transfIdx << endl;
            
            Mat img;
            if (data.useImages)
            {
                cout << data.imageNameVec[transfIdx] << endl;
                img = imread(data.imageNameVec[transfIdx], CV_LOAD_IMAGE_COLOR);
            }
            else
            { 
                img = Mat(data.imageHeight, data.imageWidth, CV_8UC3);
                img.setTo(0);
            }
            
            
            cout << transfVec[transfIdx] << endl;
            cout << "standard deviation : " << sigma << endl;
//            cout << "6 sigma : " << 6 * sigma << endl;
            
            drawPoints(img, data.detectedCornersVec[transfIdx]);
            //projected
            for (int i = 0; i < inlierVec.size(); i++)
            {
                circle(img, Point(inlierVec[i][0], inlierVec[i][1]), 9, Scalar(0, 255, 0), 5);
                
//                circle(img, Point(inlierVec[i][0], inlierVec[i][1]), 15, Scalar(0, 127, 255), 8);
            }
            for (int i = 0; i < outlierVec.size(); i++)
            {
                circle(img, Point(outlierVec[i][0], outlierVec[i][1]), 9, Scalar(0, 127, 255), 5);

//                circle(img, Point(outlierVec[i][0], outlierVec[i][1]), 15, Scalar(0, 127, 255), 8);
                
                cout << outlierVec[i].transpose() << "   err : " << errVec[i] << endl;
            }
            //detected
            
//            for (auto & pt : data.detectedCornersVec[transfIdx])
//            {
//                cross(img, pt[0], pt[1], 5, Scalar(255, 127, 0));
//            }
            
            if (data.saveOutlierImages) 
            {
                imwrite( "outliers_" + to_string(transfIdx) + ".png", img);
            }
            if (data.showOutliers) 
            {
                imshow("outliers", img);
                waitKey();
            }
        }
    }
    
    //release resources
    delete cam;
    residualFile.close();
}

