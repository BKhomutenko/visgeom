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

#include <glog/logging.h>

#include "calibration/calib_cost_functions.h"
#include "calibration/corner_detector.h"
#include "projection/generic_camera.h"
#include "projection/eucm.h"
    
bool GenericCameraCalibration::compute()
{
    //run the solver
    Solver::Options options;
//        options.check_gradients = true;
    options.gradient_check_relative_precision = 1e-5;
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
        assert(not info.constant or info.prior);
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
        if (cameraType == "eucm")
        {
            cout << "Model : EUCM" << endl;
            assert(intrinsicVec.size() == 6);
            cameraMap[name] = new EnhancedCamera(intrinsicVec.data());
        }
        else
        {
            cout << "ERROR : invalid camera model name" << endl;
            assert(false);
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
    cout << endl;
}

void GenericCameraCalibration::initGrid(ImageData & data, const ptree & node)
{   
    data.Nx = node.get<int>("object.cols");
    data.Ny = node.get<int>("object.rows");
    double sqSize = node.get<double>("object.size");
    data.board.clear();
    data.board.reserve(data.Nx * data.Ny);
    for (int i = 0; i < data.Ny; i++)
    {
        for (int j = 0; j < data.Nx; j++)
        {
           data.board.emplace_back(sqSize * j, sqSize * i, 0); 
        }
    }
    //fill up detectedCornersVec which stores all the extracted grids
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
                    ptrVec[0], intrinsicMap[data.cameraName].data());
            break;
        case 2:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], intrinsicMap[data.cameraName].data());
            break;
        case 3:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2], intrinsicMap[data.cameraName].data());
            break;
        case 4:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], intrinsicMap[data.cameraName].data());
            break;
        case 5:
            problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], ptrVec[4], intrinsicMap[data.cameraName].data());
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
    //there is a transform to initialize
    assert(transformInfoMap.find(initName) != transformInfoMap.end());
    
    //and it belongs to the transformation chain
    auto nameIter = find(data.transNameVec.begin(), data.transNameVec.end(), initName);
    assert(nameIter != data.transNameVec.end());
    
    //it is not allowed to initialize a transformation with a prior
    if (not (transformInfoMap[initName].prior or transformInfoMap[initName].initialized))
    {
        transformInfoMap[initName].initialized = true;
        
        //make sure that the data is not initializad
        if (not transformInfoMap[initName].global)
        {
            assert(sequenceTransformMap[initName].size() == 0);
        }
        
        //do the initialization
        if (not transformInfoMap[initName].global)
        {
            for (int transfIdx = 0; transfIdx < data.detectedCornersVec.size(); transfIdx++)
            {
                if (data.detectedCornersVec[transfIdx].empty())
                {
                    cout << "WARNING : " << initName << " " << transfIdx
                         << " is not initialized, no board extracted" << endl;
                    sequenceTransformMap[initName].push_back(Array6d{0, 0, 1, 0, 0, 0});    
                }
                else
                {
                    auto xi = estimateInitialGrid(data, transfIdx);
                    xi = getInitTransform(xi, initName, data, transfIdx);
                    sequenceTransformMap[initName].push_back(xi.toArray());
                }
            }
        }
        else
        {
            int transfIdx = data.getFirstExtractedIdx();
            auto xi = estimateInitialGrid(data, transfIdx);
            xi = getInitTransform(xi, initName, data, transfIdx);
            xi.toArray(globalTransformMap[initName].data());
            
            //FIXME to visualize he detected corners
            /*
            for (int i = 0; i < detectedCornersVec.size(); i++)
            {
                if (detectedCornersVec[i].empty())
                {
                    continue;
                }
                Mat8u img(800, 600);
                img.setTo(0);
                
                for (int j = 0; j < detectedCornersVec[i].size(); j++)
                {
                    img(detectedCornersVec[i][j][1], detectedCornersVec[i][j][0]) = 255;
                }
                
                
                imshow(initName, img);
                waitKey();
            }
            */
            if (data.detectedCornersVec.size() > 1) initGlobalTransform(data, initName);
        }
    }
    
    //to make sure that all the other transformations are initialized
    for (auto & x : data.transNameVec)
    {
        assert(transformInfoMap[x].prior xor transformInfoMap[x].initialized);
    }
}

void GenericCameraCalibration::addGridResidualBlocks(const ImageData & data)
{
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
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1), intrinsicPtr);
            break;
        case 1:
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], intrinsicPtr);
            break;
        case 2:
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], intrinsicPtr);
            break;
        case 3:
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2], intrinsicPtr);
            break;
        case 4:
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], intrinsicPtr);
            break;
        case 5:
            globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                    ptrVec[0], ptrVec[1], ptrVec[2],
                    ptrVec[3], ptrVec[4], intrinsicPtr);
            break;
        default:
            assert(false);
        }
        
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
        else if (dataType == "odometry")
        {
            const string transformName = dataInfo.second.get<string>("transform");
            assert(transformInfoMap.find(transformName) != transformInfoMap.end());
            assert(not transformInfoMap[transformName].global); //works only fo sequences
            
            const double errV = dataInfo.second.get<double>("err_v"); //relative error in speed
            const double errW = dataInfo.second.get<double>("err_w"); //relative error in rotation
            const double lambda = dataInfo.second.get<double>("lambda"); //relative error in rotation
            //read out the transformations
            vector<Transf> odometryVec;
            for (auto & odomItem : dataInfo.second.get_child("value"))
            {
                odometryVec.emplace_back(readTransform(odomItem.second));
            }
            
            //use the odometry as initial values
            if (dataInfo.second.get<bool>("init"))
            {
                cout << transformName << endl;
                for (auto & xx : sequenceTransformMap[transformName])
                {
                    for (auto & x : xx)
                    {
                        cout << x << "   ";
                    }
                    cout << endl;
                }
                assert(sequenceTransformMap[transformName].size() == 0);
                transformInfoMap[transformName].initialized = true;
                for (auto & xi : odometryVec)
                {
                    sequenceTransformMap[transformName].emplace_back(xi.toArray());
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
            
            assert(globalTransformMap.find(transformName) != globalTransformMap.end());
            assert(transformInfoMap[transformName].prior);
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

void GenericCameraCalibration::extractGridProjections(ImageData & data)
{
    Size patternSize(data.Nx, data.Ny);
    for (auto & fileName : data.imageNameVec)
    {
        cout << "." << flush;
        data.detectedCornersVec.emplace_back();
        Mat8u frame = imread(fileName, 0);
        if (frame.empty())
        {
            cout << fileName << " : ERROR, file not found" << endl;
            continue;
        }
        vector<cv::Point2f> centers;
        bool patternIsFound = findChessboardCorners(frame, patternSize, centers, CV_CALIB_CB_ADAPTIVE_THRESH);
        if (not patternIsFound)
        {
            cout << fileName << " : ERROR, pattern not found" << endl;
            continue;
        }
        
        if (data.checkExtraction)
        {
            drawChessboardCorners(frame, patternSize, Mat(centers), patternIsFound);
            imshow("corners", frame);
            char key = waitKey();
            if (key == 'n' or key == 'N')
            {
                cout << fileName << " : ERROR, pattern not accepted" << endl;
                continue;
            }
        }
        
        auto & cornerVec = data.detectedCornersVec.back();
        cornerVec.reserve(centers.size());
        for (auto pt : centers)
        {
            cornerVec.emplace_back(pt.x, pt.y);
        }
        
        if (data.improveDetection)
        {
            double minDist = findMinDistance(cornerVec, data.Ny, data.Nx);
            CornerDetector detector(frame, min(minDist / 2., 15.));
            detector.improveCorners(cornerVec);
        }
    }
    cout << endl;
}

Transf GenericCameraCalibration::estimateInitialGrid(const ImageData & data, const int gridIdx)
{
    auto & cornerVec = data.detectedCornersVec[gridIdx];
    
    Problem problem;
    GenericProjectionJac * costFunction = new GenericProjectionJac(cornerVec, data.board,
            cameraMap[data.cameraName], {TRANSFORM_DIRECT});
    array<double, 6> xiArr{0, 0, 1, 0, 0, 0};
    
    Vector2d v = cornerVec[1] - cornerVec[0];
    xiArr[5] = atan2(v[1], v[0]);
    
    problem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
            xiArr.data(), intrinsicMap[data.cameraName].data());
    problem.SetParameterBlockConstant(intrinsicMap[data.cameraName].data());
    Solver::Options options;
    options.max_num_iterations = 500;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    return Transf(xiArr.data());
}

void GenericCameraCalibration::computeTransforms(const ImageData & data, vector<Transf> & transfVec) const
{
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
            }
            else if (data.transStatusVec[i] == TRANSFORM_INVERSE)
            {
                xi = xi.composeInverse(getTransform(name, transfIdx));
            }
        }
    }
}

//TODO put elsewhere
void cross(Mat& img, Point pt, int size, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
{
    line(img, Point(pt.x - size, pt.y - size),
            Point(pt.x + size, pt.y + size),
            color, thickness, lineType, shift);
    line(img, Point(pt.x - size, pt.y + size),
            Point(pt.x + size, pt.y - size),
            color, thickness, lineType, shift);
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
        
        double stdAcc;
        Vector2dVec inlierVec, outlierVec;
        vector<double> errVec;
        for (int i = 0; i < projectedVec.size(); i++)
        {
            Vector2d err = data.detectedCornersVec[transfIdx][i] - projectedVec[i];
            residualFile << err.transpose() << "   " << projectedVec[i].transpose() << endl;
            stdAcc += err.squaredNorm();
        }
        
        double sigma = sqrt(stdAcc / (projectedVec.size() - 1));
        for (int i = 0; i < projectedVec.size(); i++)
        {
            Vector2d err = data.detectedCornersVec[transfIdx][i] - projectedVec[i];
            double errNorm = err.norm();
            if (errNorm < 3 * sigma and errNorm < 1) // 1px is fo the case when all the points are off
            {                                        // we are looking for the sub-pixel precision
                inlierVec.push_back(projectedVec[i]);
            }
            else
            {
                outlierVec.push_back(projectedVec[i]);
                errVec.push_back(err.norm());
            }
        }
        
        if (data.showOutliers and not outlierVec.empty()) //TODO make a flag to display the outliers
        {
            cout << data.imageNameVec[transfIdx] << endl;
            cout << "standard deviation : " << sigma << endl;
            cout << "3 sigma : " << 3 * sigma << endl;
            Mat img = imread(data.imageNameVec[transfIdx], CV_LOAD_IMAGE_COLOR);
            
            //projected
            for (int i = 0; i < inlierVec.size(); i++)
            {
                circle(img, Point(inlierVec[i][0], inlierVec[i][1]), 8, Scalar(0, 255, 0), 3);
            }
            for (int i = 0; i < outlierVec.size(); i++)
            {
                circle(img, Point(outlierVec[i][0], outlierVec[i][1]), 8, Scalar(0, 127, 255), 3);
                cout << outlierVec[i] << "   err : " << errVec[i] << endl;
            }
            
            //detected
            for (auto & pt : data.detectedCornersVec[transfIdx])
            {
                cross(img, Point(pt[0], pt[1]), 5, Scalar(255, 127, 0));
            }
            
            imshow("outliers", img);
            waitKey();
        }
    }
    
    //release resources
    delete cam;
    residualFile.close();
}

