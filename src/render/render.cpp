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

#include "render/render.h"

#include "timer.h"

#include "projection/eucm.h"

#include "render/background.h"
#include "render/plane.h"

Renderer::Renderer(const ptree & params) :
    _width(params.get<int>("width")),
    _height(params.get<int>("height")),
    _camera(NULL)
{
    Size size(_width, _height);
    _idxMat.create(size);
    _uMat.create(size);
    _vMat.create(size);
    _depthMat.create(size);
    
    _objectVec.push_back( new Background(params.get_child("background")) );
    
    for (auto & objParams : params.get_child("plains"))
    {
        _objectVec.push_back( new Plane(objParams.second) );
    }
    
    //TODO init the objects
}
    
Renderer::~Renderer()
{
    if (_camera != NULL)
    {
        delete _camera;
    }
    for (auto obj : _objectVec)
    {
        delete obj;
    }
}

/*
    Algorithm:
    -fill up the buffers using object->intersect
    -for each pixel compute local basis
    -sample textures
*/    
void Renderer::setCameraTransform(const Transf & xi) 
{ 
    _xiCam = xi; 
}

void Renderer::setCamera(const ICamera * camera) 
{
    if (_camera != NULL)
    {
        delete _camera;
    }
    _camera = camera->clone();
}

void Renderer::fillBuffers() 
{
    _idxMat.setTo(-1);
    _depthMat.setTo(1e6);
    Matrix3d R = _xiCam.rotMat();
    for (int v = 0; v < _height; v++)
    {
        for (int u = 0; u < _width; u++)
        {
            Vector3d dir;
            if (not  _camera->reconstructPoint(Vector2d(u, v), dir) ) continue;
            dir = R * dir;
            for (int idx = 0; idx < _objectVec.size(); idx++)
            {
                Vector2d uv; //texture coordinates
                double depth;
                if (not _objectVec[idx]->intersection(_xiCam.trans(), dir, uv, depth)) continue;
                
                if (_depthMat(v, u) > depth)
                {
                    _depthMat(v, u) = depth;
                    _idxMat(v, u) = idx;
                    _uMat(v, u) = uv[0];
                    _vMat(v, u) = uv[1];
                }
            }
        }
    }
}

void Renderer::fillImage(Mat8u & dst) 
{
    Timer timer;
    
    vector<double> tVec(_objectVec.size(), 0);
    vector<int> pxCount(_objectVec.size(), 0);
    dst.create(_height, _width);
    for (int v = 0; v < _height; v++)
    {
        for (int u = 0; u < _width; u++)
        {
            int idx = _idxMat(v, u);
            if (idx == -1) continue;
            
            Vector2d pt(_uMat(v, u), _vMat(v, u));
            
            //first, find a local basis
            Matrix2d basis;
            basis << 0, 0, 0, 0;
            
            //along u
            int uCount = 0;
            double uMin = pt[0], uMax = pt[0];
            double vMin = pt[1], vMax = pt[1];
            if (u > 0 and _idxMat(v, u - 1) == idx)
            {
                uCount++;
                uMin = _uMat(v, u - 1);
                vMin = _vMat(v, u - 1);
            }
            if (u < _width - 1 and _idxMat(v, u + 1) == idx)
            {
                uCount++;
                uMax = _uMat(v, u + 1);
                vMax = _vMat(v, u + 1);
            }
            
            if (uCount != 0)
            {
                basis(0, 0) = (uMax - uMin) / uCount;
                basis(1, 0) = (vMax - vMin) / uCount;
            }
            
            //along v
            int vCount = 0;
            uMin = pt[0], uMax = pt[0];
            vMin = pt[1], vMax = pt[1];
            if (v > 0 and _idxMat(v - 1, u) == idx)
            {
                vCount++;
                uMin = _uMat(v - 1, u);
                vMin = _vMat(v - 1, u);
            }
            if (v < _height - 1 and _idxMat(v + 1, u) == idx)
            {
                vCount++;
                uMax = _uMat(v + 1, u);
                vMax = _vMat(v + 1, u);
            }
            
            if (vCount != 0)
            {
                basis(0, 1) = (uMax - uMin) / vCount;
                basis(1, 1) = (vMax - vMin) / vCount;
            }
            
            timer.reset();            
            dst(v, u) = _objectVec[idx]->sample(pt, basis);
            tVec[idx] += timer.elapsed();
            pxCount[idx]++;
        }
    }
    if (false) //FIXME make verbose
    for (int i = 0; i < tVec.size(); i++)
    {
        cout << setw(15) << pxCount[i] << setw(15) << tVec[i] << setw(15) << tVec[i] / pxCount[i] << endl;
    }
}


