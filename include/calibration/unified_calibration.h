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
#include "utils/json.h"

//TODO make a separate header


#include <glog/logging.h>

#include "calib_cost_functions.h"



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
    
    void initGlobalTransform(const string & name);
    
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

/* TODO add the residual analysis
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
    }
*/

