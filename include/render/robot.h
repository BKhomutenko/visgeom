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
#include "eigen.h"
#include "geometry/geometry.h"
#include "projection/generic_camera.h"

#include "render/render.h"

class VirtualRobot
{
public:
    VirtualRobot(const ptree & params);
        
    virtual ~VirtualRobot();
  
    
    void setVelocity(const Vector3d v, const Vector3d omega); //FIXME implement a separate data structure

    void simulationStep(const double timeStep);
    
    void getImage(Mat8u & dst, int cameraIdx = 0);
    
    void fillBuffers();
    
    void fillImage(Mat8u & dst);
    
    vector<ICamera*> _cameraVec;
    vector<Transf> _xiBaseCamVec;
    Transf _xiOrigBase;
    
    Renderer _renderDevice;
};

