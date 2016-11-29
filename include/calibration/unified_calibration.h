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
    bool constant;
    bool initialized;  
};

//TODO choose the model in the parameter file, not as atemplate parameter
class GenericCameraCalibration
{
protected:
    map<string, vector<double>> intrinsicMap;
    map<string, Array6d> globalTransformMap;
    map<string, vector<Array6d>> sequenceTransformMap;
    map<string, TransformInfo> transformInfoMap;
    map<string, ICamera*> cameraMap;
    Problem globalProblem;
    
    //temporary variables
    int Nx, Ny;
    ptree root;
    vector<Vector2dVec> gridExtractionVec;
    vector<string> transNameVec;
    vector<TransformationStatus> transStatusVec;
    Vector3dVec grid;
    string cameraName;
    
    void parseTransforms();
    
    void parseCameras();
    
    void parseData();
    
    void initTransformChainInfo(const ptree & node);
    
    void initGrid(const ptree & node);
    
    void initTransforms(const ptree & node);
    
    void getInitTransform(Transformation<double> & xi, const string & initName, int gridIdx = 0);
    
    void addGridResidualBlocks();
    
public:
    virtual ~GenericCameraCalibration() 
    { 
        for (auto & x : cameraMap)
        {
            delete x.second;
        }
    }
    
    bool compute();
    
    //returns the first transform if it is a sequence
    Array6d & getTransformData(const string & name, int idx = 0)
    {
        if (transformInfoMap[name].global)  return globalTransformMap[name];
        else return sequenceTransformMap[name][idx];
    }
    
    Transformation<double> getTransform(const string & name, int idx = 0)
    {
        return Transformation<double>(getTransformData(name, idx).data());
    }
    
    
    
    bool addResiduals(const string & infoFileName)
    {
        read_json(infoFileName, root);
        parseTransforms();
        parseCameras();
        parseData();
    }
    
    bool extractGridProjection(const string & fileName, Vector2dVec & projection, bool checkExtraction);

    Transformation<double> estimateInitialGrid(const string & cameraName,
            const Vector2dVec & projection, const Vector3dVec & grid);

};

