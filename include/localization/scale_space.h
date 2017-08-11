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

/*
Scale space for multiscale optimization
*/

#pragma once

#include "std.h"
#include "ocv.h"
#include "io.h"

class BinaryScalSpace
{
public:
    BinaryScalSpace(int numScales = 1, bool withGradient = false) : 
            activeScaleIdx(0), 
            gradientOn(withGradient)
    {
        scale = 1;
        assert(numScales > 0);
        imgVec.resize(numScales);
        if (gradientOn) resizeGradient();
    }
    
    int uConv(int x) const { return x * scale; }
    int vConv(int y) const { return y * scale; }
    
    void setGradient(bool val)
    {
        gradientOn = val; 
        if (gradientOn) resizeGradient();
    }
    
    void setNumberScales(int val)
    {
        assert(val > 0);
        imgVec.resize(val);
        if (gradientOn) resizeGradient();
    }
    
    void generate(const Mat8u & img)
    {
        img.convertTo(imgVec[0], CV_32F);
        propagate();
    }

    void generate(const Mat32f & img)
    {
        img.copyTo(imgVec[0]);
        propagate();
    }

    const Mat32f & get() const { return imgVec[activeScaleIdx]; }
    const Mat32f & getGradU() const { return gradUVec[activeScaleIdx]; }
    const Mat32f & getGradV() const { return gradVVec[activeScaleIdx]; }
    
    int size() const { return imgVec.size(); }
    
    int scaleByIdx(int idx) const { return (1 << idx); }
    
    void setActiveScale(int idx) 
    { 
        scale = scaleByIdx(idx);
        activeScaleIdx = idx;
    }
    
    int getActiveScale() 
    { 
        return scale;
    }
    
    int getActiveIdx()
    {
        return activeScaleIdx;
    }
    
private:

    void resizeGradient()
    {
        gradUVec.resize(size());
        gradVVec.resize(size());
    }
    
    void propagate()
    {
        if (gradientOn)
        {
            Sobel(imgVec[0], gradUVec[0], CV_32F, 1, 0, 3, 1./8);
            Sobel(imgVec[0], gradVVec[0], CV_32F, 0, 1, 3, 1./8);
        }
        for (int i = 1; i < imgVec.size(); i++)
        {
            /*Mat32f blured;
            GaussianBlur(imgVec[i - 1], blured, Size(3, 3), 0, 0);
            imgVec[i].create(Size(blured.cols/2, blured.rows/2));
            for (int v = 0; v < imgVec[i].rows; v++)
            {
                for (int u = 0; u < imgVec[i].cols; u++)
                {
                    imgVec[i](v, u) = blured(2*v, 2*u);
                }
            }*/
            
            imgVec[i].create(Size(imgVec[i - 1].cols/2, imgVec[i - 1].rows/2));
            imgVec[i].setTo(0);
            for (int v = 0; v < imgVec[i - 1].rows; v++)
            {
                int vs = min((int)round(v / 2.), imgVec[i].rows - 1);
                for (int u = 0; u < imgVec[i - 1].cols; u++)
                {
                    int us = min((int)round(u / 2), imgVec[i].cols - 1);
                    imgVec[i](vs, us) += imgVec[i - 1](v, u);
                }
            }
            imgVec[i] *= 0.25;
            if (gradientOn)
            {
                Sobel(imgVec[i], gradUVec[i], CV_32F, 1, 0, 3, 1./8);
                Sobel(imgVec[i], gradVVec[i], CV_32F, 0, 1, 3, 1./8);
            }
        }
    }
    
    std::vector<Mat32f> imgVec;
    std::vector<Mat32f> gradUVec;
    std::vector<Mat32f> gradVVec;
    int scale;
    int activeScaleIdx;
    bool gradientOn;
};
