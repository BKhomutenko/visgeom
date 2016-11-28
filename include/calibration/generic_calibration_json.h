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
        options.gradient_check_relative_precision = 1e-2;
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
                auto nameIter = find(transNameVec.begin(), transNameVec.end(), initName);
                assert(nameIter != transNameVec.end());
                assert(transformInfoMap.find(initName) != transformInfoMap.end());
                
                //it is not allowed to initialize a transformation with a prior
                assert(transformInfoMap[initName].prior == false);
                assert(transformInfoMap[initName].initialized == false);
                if (not transformInfoMap[initName].global)
                {
                    assert(sequenceTransformMap[initName].size() == 0);
                }
                initialize = true;
                transformInfoMap[initName].initialized = true;
            }
            
            //to make sure that all the other transformations are initialized
            for (auto & x : transNameVec)
            {
                assert(transformInfoMap[x].prior xor
                     transformInfoMap[x].initialized);
            }
            
            //process images, init the transform
            //TODO separate globalProblem initialization and grid extraction:
            //1 - create a grid vector which stores all the extracted grids
            //2 - initialize the transforms and 
            //3 - initialize globalProblem
            Vector2dVec projection;
            int sequenceIdx = -1;
            for (auto & x : dataInfo.second.get_child("images.names"))
            {
                string filename = x.second.get_value<string>();
                cout << ".";
                sequenceIdx++;
                if (not extractGridProjection(prefix + filename, projection, checkExtraction))
                {
                    cout << "WARNING : GRID NOT EXTRACTED" << endl;
                    if (initialize and not transformInfoMap[initName].global) 
                    {   
                        //to make sure that the number of transforms corresponds to the number of images
                        // even though some of them are not initialized
                        sequenceTransformMap[initName].push_back({0, 0, 1, 0, 0, 0});
                    }
                    continue;
                }
                
                //TODO properly initialize a global transformation
                // make two separate functions
                // initGlobal, initSequence
                if (initialize)
                {
                    auto xi = estimateInitialGrid(cameraMap[cameraName], projection, grid);
                    if (transformInfoMap[initName].global) cout << xi << endl;
                    auto constTransIter = constTransVec.begin();
                    for (int i = 0; i < transNameVec.size(); i++)
                    {
                        const string & name = transNameVec[i];
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
                        if (name == initName)
                        {
                            if (transStatusVec[i] == TRANSFORM_INVERSE)
                            {
                                xi = xi.inverse();
                            }
                            break;
                        }
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
                        if (transformInfoMap[initName].global)
                        {
                            cout << xi << endl;
                            cout << getTransform(name, sequenceIdx) << endl;
                        }
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
                        auto & arr = getTransformData(name, sequenceIdx);
                        ptrVec.push_back(arr.data());
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
            cout << endl; // to close the sequence of "." TODO make this output properly
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
        options.max_num_iterations = 500;
        Solver::Summary summary;
        Solve(options, &problem, &summary);
        return Transformation<double>(xi.data());
    }

};

