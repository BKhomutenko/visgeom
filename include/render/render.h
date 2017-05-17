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

#include "json.h"
#include "std.h"
#include "ocv.h"
#include "eigen.h"
#include "geometry/geometry.h"
#include "projection/generic_camera.h"

#include "render/object.h"

class Renderer
{
public:
    Renderer(const ptree & params);
        
    virtual ~Renderer();
  
    void setCameraTransform(const Transf & xi);
    void setCamera(const ICamera * camera);
    
    //TODO put to a cpp file
    void fillBuffers();
    
    void fillImage(Mat8u & dst);
    
    vector<IObject * > _objectVec;
    ICamera * _camera;
    
    //TODO init all the matrices
    Mat16s _idxMat;
    Mat32f _uMat, _vMat, _depthMat;
    Transf _xiCam;
    int _width, _height;
};

