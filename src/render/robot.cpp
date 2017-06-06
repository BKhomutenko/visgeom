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


#include "render/robot.h"

VirtualRobot::VirtualRobot(const ptree & params)
{
    for (auto & element : params.get_child("cameras"))
    {
        auto camParams = element.second;
        vector<double> intrinsics = readVector(camParams.get_child("intrinsics"));
        _cameraVec.push_back( new EnhancedCamera(intrinsics.data()) );
        _xiBaseCamVec.push_back( readTransform(camParams.get_child("extrinsics")) );
    }
}
        
VirtualRobot::~VirtualRobot()
{
    for (auto & ptr : _cameraVec)
    {
        delete ptr;
    }
}

void VirtualRobot::setVelocity(const Vector3d v, const Vector3d omega)
{
    _v = v;
    _omega = omega;
}

void VirtualRobot::simulationStep(const double timeStep)
{
    double omegaAbs = _omega.norm();
    if (omegaAbs * timeStep < 1e-2) //TODO revise the literal
    {
        //half rotation
        Matrix3d R2 = rotationMatrix(_omega * 0.5 * timeStep);
        Vector3d vEffective = R2 * _v;
        
        Transf zeta(vEffective * timeStep, _omega * timeStep);
        _xiOrigBase = _xiOrigBase.compose(zeta);
    }
    else
    {
        //rotation radius vector
        Vector3d r = _v.cross(_omega) / (omegaAbs * omegaAbs); 
        
        Vector3d vAxial = _v.dot(_omega) / omegaAbs;
        Matrix3d R = rotationMatrix(_omega * timeStep);
        
        Vector3d trans = vAxial * timeStep - r + R * r;
        
        Transf zeta(trans, _omega * timeStep);
        _xiOrigBase = _xiOrigBase.compose(zeta);
    }
}

void VirtualRobot::renderImage(RenderDevice & device, Mat8u & dst, int cameraIdx = 0)
{
    if (cameraIdx >= _cameraVec.size()) throw std::runtime_error("Robot : camera index is out of range");
    device.setCameraTransformation(_xiOrigBase.compose(_xiBaseCamVec[cameraIdx]));
    device.setCamera(_cameraVec[idx]);
    device.render(dst);
}

void VirtualRobot::fillBuffers();

void VirtualRobot::fillImage(Mat8u & dst);



