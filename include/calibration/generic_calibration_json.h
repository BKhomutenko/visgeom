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

struct TransformInfo {
    bool global;
    bool prior;
    bool initialized;  
};

//TODO choose the model in the parameter file, not as atemplate parameter
class GenericCameraCalibration
{
protected:
    map<string, vector<double>> intrinsicMap;
    map<string, Array6d> globalTransformMap;
    map<string, list<Array6d>> sequenceTransformMap;
    map<string, TransformInfo> transformInfoMap;
    map<string, ICamera*> cameraMap;
    Problem globalProblem;
    
    //temporary variables
    int Nx, Ny;
    ptree root;
    
public:
    //TODO implement access to the results
    
    virtual ~GenericCameraCalibration()
    { }
    
    bool compute()
    {
        //run the solver
        Solver::Options options;
//        options.check_gradients = true;
        options.max_num_iterations = 250;
        options.function_tolerance = 1e-10;
        options.gradient_tolerance = 1e-10;
        options.parameter_tolerance = 1e-10;
        options.logging_type = ceres::SILENT;
//        options.minimizer_progress_to_stdout = true;
        Solver::Summary summary;
        Solve(options, &globalProblem, &summary);
        cout << summary.FullReport() << endl;
        
        cout << "Intrinsic parameters :" << endl;
        for (auto & x : intrinsicMap)
        {
            cout << x.first << " : ";
            for (int i = 0; i < cameraMap[x.first]->numParams(); i++) //FIXME store the intrinsic size
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
                cout << i++ << " : " << Transformation<double>(x.data()) << endl;
            }
        }
        
        cout << "Global extrinsic parameters :" << endl;
        for (auto & x : globalTransformMap)
        {
            cout << x.first << " : " << Transformation<double>(x.second.data()) << endl;
        }
        
    }
    
    void parseTransforms()
    {
        for (auto & transInfo : root.get_child("transformations"))
        {
            string name = transInfo.second.get<string>("name");
            transformInfoMap[name] = TransformInfo();
            auto & info = transformInfoMap[name];
            info.global = transInfo.second.get<bool>("global");
            info.prior = transInfo.second.get<bool>("prior");
            info.initialized = false;
            
            if (info.global)
            {
                globalTransformMap[name] = Array6d();
                if (info.prior)
                {
                    int i = 0;
                    for (auto & x : transInfo.second.get_child("value"))
                    {
                        globalTransformMap[name][i++] = x.second.get_value<double>();
                    }
                }
            }
            else 
            {
                sequenceTransformMap[name] = list<Array6d>();
                auto & valVec = sequenceTransformMap[name];
                if (info.prior)
                {
                    for (auto & val : transInfo.second.get_child("value"))
                    {
                        valVec.emplace_back();
                        int i = 0;
                        for (auto & x : val.second)
                        {
                            valVec.back()[i++] = x.second.get_value<double>();
                        }
                    }
                }
            }
        }
    }
    
    void parseCameras()
    {
        for (auto & cameraInfo : root.get_child("cameras"))
        {
            string name = cameraInfo.second.get<string>("name");
            intrinsicMap[name] = vector<double>();
            auto & intrinsicVec = intrinsicMap[name];
            string cameraType = cameraInfo.second.get<string>("type");
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
    
    //returns the first transform if it is a sequence
    Array6d & getTransformData(const string & name, int idx = 0)
    {
        if (transformInfoMap[name].global)  return globalTransformMap[name];
        else
        {
            auto iter = sequenceTransformMap[name].begin();
            advance(iter, idx);
            return *iter;
        }
    }
    
    Transformation<double> getTransform(const string & name, int idx = 0)
    {
        return Transformation<double>(getTransformData(name, idx).data());
    }
    
    void parseData()
    {
        for (auto & dataInfo : root.get_child("data"))
        {
            //TODO add different init for a different type
            vector<string> transNameVec;
            vector<TransformationStatus> transStatusVec;
            vector<Transformation<double>> constTransVec;
            for (auto & transInfo : dataInfo.second.get_child("transform_chain"))
            {
                transNameVec.push_back(transInfo.second.get<string>("name"));
                
                if (transformInfoMap.find(transNameVec.back()) == transformInfoMap.end())
                {
                    //the transformation is constant
                    transStatusVec.push_back(TRANSFORM_CONSTANT);
                    vector<double> xiData;
                    xiData.reserve(6);
                    for (auto & x : transInfo.second.get_child("value"))
                    {
                        xiData.push_back(x.second.get_value<double>());
                    }
                    assert(xiData.size() == 6);
                    constTransVec.emplace_back(xiData.data());
                }
                else
                {
                    if (transInfo.second.get<bool>("direct")) 
                    {
                        transStatusVec.push_back(TRANSFORM_DIRECT);
                    }
                    else transStatusVec.push_back(TRANSFORM_INVERSE);
                }
            }
            
            string cameraName = dataInfo.second.get<string>("camera");
            
            Nx = dataInfo.second.get<int>("object.cols");
            Ny = dataInfo.second.get<int>("object.rows");
            double sqSize= dataInfo.second.get<double>("object.size");
            Vector3dVec grid;
            grid.reserve(Nx*Ny);
            for (int i = 0; i < Ny; i++)
            {
                for (int j = 0; j < Nx; j++)
                {
                   grid.emplace_back(sqSize * j, sqSize * i, 0); 
                }
            }
            
            bool checkExtraction = dataInfo.second.get<bool>("parameters.check_extraction");
            string initName = dataInfo.second.get<string>("init");
            string prefix = dataInfo.second.get<string>("images.prefix");
            
            bool initialize = false;
            if (initName != "none")
            {
                //there is a transform to initialize
                assert(transformInfoMap.find(initName) != transformInfoMap.end());
                assert(transformInfoMap[initName].prior == false);
                assert(transformInfoMap[initName].initialized == false);
                if (not transformInfoMap[initName].global)
                {
                    assert(sequenceTransformMap[initName].size() == 0);
                }
                initialize = true;
                transformInfoMap[initName].initialized = true;
            }
            
            for (auto & x : transNameVec)
            {
//                cout << x << endl;
//                cout << transformInfoMap[x].prior << " " <<
//                     transformInfoMap[x].initialized << endl;
                assert(transformInfoMap[x].prior xor
                     transformInfoMap[x].initialized);
            }
            
            //process images, init the transform
            Vector2dVec projection;
            int sequenceIdx = -1;
            for (auto & x : dataInfo.second.get_child("images.names"))
            {
                //TODO make possible const xiImage
                string filename = x.second.get_value<string>();
                cout << filename << endl;
                sequenceIdx++;
                if (not extractGridProjection(prefix + filename, projection, checkExtraction))
                {
                    if (not transformInfoMap[initName].global)
                    {
                        sequenceTransformMap[initName].push_back({0, 0, 1, 0, 0, 0});
                    }
                    continue;
                }
                
                
                if (initialize)
                {
    //                cout << "    Estimate initial extrinsic" << endl;
                    auto xi = estimateInitialGrid(cameraMap[cameraName], projection, grid);
                    auto constTransIter = constTransVec.begin();
                    for (int i = 0; i < transNameVec.size(); i++)
                    {
                        const string & name = transNameVec[i];
    //                    cout << "        " << name << endl;
                        if (name == initName) break;
                        else if (transStatusVec[i] == TRANSFORM_CONSTANT)
                        {
                            xi = (*constTransIter).inverseCompose(xi);
                            constTransIter++;
                        }
                        else if (transStatusVec[i] == TRANSFORM_DIRECT)
                        {
                            xi = getTransform(name, sequenceIdx).inverseCompose(xi);
                        }
                        else if (transStatusVec[i] == TRANSFORM_INVERSE)
                        {
                            xi = getTransform(name, sequenceIdx).compose(xi);
                        }
                    }
                    constTransIter = constTransVec.end();
                    for (int i = transNameVec.size() - 1; i >= 0; i--)
                    {
                        const string & name = transNameVec[i];
    //                    cout << "        " << name << endl;
                        if (name == initName) break;
                        else if (transStatusVec[i] == TRANSFORM_CONSTANT)
                        {
                            constTransIter--;
                            xi = xi.composeInverse(*constTransIter);
                        }
                        else if (transStatusVec[i] == TRANSFORM_DIRECT)
                        {
                            xi = xi.composeInverse(getTransform(name, sequenceIdx));
                        }
                        else if (transStatusVec[i] == TRANSFORM_INVERSE)
                        {
                            xi = xi.compose(getTransform(name, sequenceIdx));
                        }
                    }
                    
                    if (transStatusVec.back() == TRANSFORM_INVERSE)
                    {
                        xi = xi.inverse();
                    }
                    
                    if (transformInfoMap[initName].global)
                    {
                        xi.toArray(globalTransformMap[initName].data());
                        initialize = false;
                    }
                    else
                    {
                        sequenceTransformMap[initName].push_back(xi.toArray());
                    }
                }
                
                
                // make the vector of pointers to the transformation data
                cout << "    Initialisze pointers" << endl;
                vector<double*> ptrVec;
                for (int i = 0; i < transNameVec.size(); i++)
                {
                    const string & name = transNameVec[i];
                    if (transStatusVec[i] == TRANSFORM_CONSTANT) continue;
                    else if (transformInfoMap[name].global)
                    {
                        ptrVec.push_back(globalTransformMap[name].data());
                    }
                    else
                    {
                        ptrVec.push_back(sequenceTransformMap[name].back().data());
                    }
                }

                //add a residual
                GenericProjectionJac * costFunction = new GenericProjectionJac(projection, grid,
                                            cameraMap[cameraName], transStatusVec, constTransVec,
                                            true);
                                            
                switch (ptrVec.size())
                {
                case 0:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                         intrinsicMap[cameraName].data());
                    break;
                case 1:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], intrinsicMap[cameraName].data());
                    break;
                case 2:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], intrinsicMap[cameraName].data());
                    break;
                case 3:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2], intrinsicMap[cameraName].data());
                    break;
                case 4:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3], intrinsicMap[cameraName].data());
                    break;
                case 5:
                    globalProblem.AddResidualBlock(costFunction, new SoftLOneLoss(1),
                            ptrVec[0], ptrVec[1], ptrVec[2],
                            ptrVec[3], ptrVec[4], intrinsicMap[cameraName].data());
                    break;
                }
            }
        }
    }
    
    
    bool addResiduals(const string & infoFileName)
    {
        read_json(infoFileName, root);
        parseTransforms();
        parseCameras();
        parseData();
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

    Transformation<double> estimateInitialGrid(const ICamera * camera,
            const Vector2dVec & projection, const Vector3dVec & grid)
    {
        Problem problem;
        GenericProjectionJac * costFunction = new GenericProjectionJac(projection, grid,
                                        camera, {TRANSFORM_DIRECT}, {}, false);
        array<double, 6> xi{0, 0, 1, 0, 0, 0};
        
        Vector2d v = projection[1] - projection[0];
        xi[5] = atan2(v[1], v[0]);
        
        problem.AddResidualBlock(costFunction, new SoftLOneLoss(1), xi.data());
        
        Solver::Options options;
        options.max_num_iterations = 250;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
//        cout << summary.BriefReport() << endl;
        return Transformation<double>(xi.data());
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

