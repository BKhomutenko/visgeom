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
Cost functions for localization based on photometric data and mutual information
*/

#include "localization/cost_functions.h"

#include "std.h"
#include "eigen.h"
#include "ocv.h" //TODO replace here Mat32f by a pointer
#include "ceres.h"
#include "io.h"

#include "geometry/geometry.h"
#include "camera/generic_camera.h"
#include "utils/jacobian.h"

/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera
*/
    
bool PhotometricCostFunction::Evaluate(double const * const * parameters,
        double * residual, double ** jacobian) const
{
    Transformation<double> T12(parameters[0]);
    
    // point cloud in frame 2
    vector<Vector3d> transformedPoints;
    T12.inverseTransform(_dataPack.cloud, transformedPoints);
    
    // init the image interpolation
    ceres::BiCubicInterpolator<Grid2D> imageInterpolator(_imageGrid);
    
    bool computeJac = (jacobian != NULL and jacobian[0] != NULL);
    
    if (computeJac)
    {
        // L_uTheta
        CameraJacobian jacobianCalculator(_camera, T12);
        for (int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2d pt;
            if (_camera->projectPoint(transformedPoints[i], pt)) 
            {
                double f;
                // image interpolation and gradient
                Covector2d grad;
                imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale,
                        &f, &grad[1], &grad[0]);
                grad /= _scale;  // normalize according to the scale
                
                jacobianCalculator.dfdxi(transformedPoints[i], grad, jacobian[0] + i*6);
                for (int j = 0; j < 6; j++) *(jacobian[0] + i*6 + j) /= (_dataPack.valVec[i] + 1); 
                
                residual[i] = (f - _dataPack.valVec[i]) / (_dataPack.valVec[i] + 1);
            }
            else
            {
                residual[i] = 0.;
                for (int j = 0; j < 6; j++)
                {
                    jacobian[0][i*6 + j] = 0.;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < transformedPoints.size(); i++)
        {
            Vector2d pt;
            if (_camera->projectPoint(transformedPoints[i], pt)) 
            {
                double f;
                imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale, &f);
                residual[i] = (f - _dataPack.valVec[i]) / (_dataPack.valVec[i] + 1);
            }
            else
            {
                residual[i] = 0.;
            }
        }
    }
    return true;
}
    
    
bool MutualInformation::Evaluate(double const * parameters,
        double * cost, double * gradient) const
{
    for (int i = 0; i < 6; i++)
    {
        if (std::isnan(parameters[i]) or std::isinf(parameters[i])) return false;
    }
    Transformation<double> T12(parameters);
    
    // point cloud in frame 2
    vector<Vector3d> transformedPoints;
    T12.inverseTransform(_dataPack.cloud, transformedPoints);
    // init the image interpolation
    ceres::BiCubicInterpolator<Grid2D> imageInterpolator(_imageGrid);
    
    bool computeGrad = (gradient != NULL);
    
    vector<double> valVec2(transformedPoints.size(), 0);
    vector<Covector2d> gradVec;
    if (computeGrad)
    {
        gradVec.resize(transformedPoints.size(), Covector2d(0, 0));
    }
    *cost = 0;
    for (int i = 0; i < transformedPoints.size(); i++)
    {
        Vector2d pt;
        auto & X = transformedPoints[i];
        if (_camera->projectPoint(transformedPoints[i], pt)) 
        {
            if (computeGrad)
            {
                double & f = valVec2[i];
                Covector2d & grad = gradVec[i];
                // image interpolation and gradient
                imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale,
                        &f, &grad[1], &grad[0]);
                grad /= _scale;  // normalize according to the scale
            }
            else 
            {
                double & f = valVec2[i];
                imageInterpolator.Evaluate(pt[1] / _scale, pt[0] / _scale, &f);
            }
            
        }
        
    }
    vector<double> hist12 = computeHist2d(_dataPack.valVec, valVec2);
    vector<double> hist2 = reduceHist(hist12);
    vector<double> log12Vec(_numBins * _numBins, 0);
    
    // compute the cost
    *cost = 0;
    for (int idx2 = 0; idx2 < _numBins; idx2++)
    {
        for (int idx1 = 0; idx1 < _numBins; idx1++)
        {
            const double & p12 = hist12[idx2 * _numBins + idx1];
            const double log12 = log(p12 / hist2[idx2] / _hist1[idx1]);
            if (p12 > 0)
            {
                log12Vec[idx2 * _numBins + idx1] = log12;
                *cost -= p12*log12;
            }
        }
    }
    // compute the gradient
    if (computeGrad)
    {
        Map<Covector6d> dMIdxi(gradient);
        dMIdxi << 0, 0, 0, 0, 0, 0;
        // L_uTheta
        CameraJacobian jacobianCalculator(_camera, T12);
        for (int i = 0; i < transformedPoints.size(); i++)
        {
            Covector6d dfdxi;
            jacobianCalculator.dfdxi(transformedPoints[i], gradVec[i], dfdxi.data());
            
            // dP/df and dMI/dP
            int idx11, idx12;
            int idx21, idx22;
            double dPdf, share1;
            computeShareDerivative(valVec2[i], idx21, idx22, dPdf);
            computeShares(_dataPack.valVec[i], idx11, idx12, share1);
            double dMIdP = 0;
           
            if (idx22 != -1)
            {
                if (idx12 != -1)
                {
                    dMIdP = log12Vec[idx21 * _numBins + idx11] * share1
                        + log12Vec[idx21 * _numBins + idx12] * (1 - share1)
                        - log12Vec[idx22 * _numBins + idx11] * share1
                        - log12Vec[idx22 * _numBins + idx12] * (1 - share1);
                }
                else
                {
                    dMIdP = log12Vec[idx21 * _numBins + idx11]
                        - log12Vec[idx22 * _numBins + idx11];
                }
            }
            double dMIdf = dMIdP * _increment * dPdf;
            dMIdxi -= dMIdf * dfdxi;
        }
    }
    return true;
}
    
void MutualInformation::computeShareDerivative(double val, int & idx1, int & idx2, double & der) const
{
    double scaledVal = val / _histStep;
    idx1 = round(scaledVal);
    
    double tail = abs(idx1 - scaledVal);
    if (idx1 < 0)
    { 
        idx1 = 0;
        idx2 = -1;
        der = 0;
    }
    else if (idx1 >= _numBins)
    { 
        idx1 = _numBins - 1;
        idx2 = -1;
        der = 0;
    }
    else
    {
        // der =  24 * (1 - 1.25 * tail) * tail * tail / _histStep;
        der =  4 * tail / _histStep;
        if (scaledVal > idx1 and idx1 < _numBins - 1)
        {
            der = -der;
            idx2 = idx1 + 1;
        }
        else if (scaledVal < idx1 and idx1 > 0)
        {
            idx2 = idx1 - 1;
        }
        else
        {
            idx2 = -1;
            der = 0;
        }
    }
}
    
void MutualInformation::computeShares(double val, int & idx1, int & idx2, double & share) const
{
    double scaledVal = val / _histStep;
    idx1 = round(scaledVal);
    double tail = abs(idx1 - scaledVal);
    if (idx1 < 0)
    { 
        idx1 = 0;
        idx2 = -1;
        share = 0;
    }
    else if (idx1 >= _numBins)
    { 
        idx1 = _numBins - 1;
        idx2 = -1;
        share = 0;
    }
    else
    {
        //  share = 1. - 8 * (tail * tail* tail) * (1 - tail);
        share = 1. - 2 * tail * tail;
        if (scaledVal > idx1 and idx1 < _numBins - 1)
        {
            idx2 = idx1 + 1;
        }
        else if (scaledVal < idx1 and idx1 > 0)
        {
            idx2 = idx1 - 1;
        }
        else idx2 = -1;
    }
}
    
vector<double> MutualInformation::computeHist(const vector<double> & valVec) const
{
    vector<double> hist(_numBins, 0);
    for (auto & x : valVec)
    {
        int idx1, idx2;
        double share;
        computeShares(x, idx1, idx2, share);
        if (idx2 != -1)
        {
            hist[idx2] += _increment * (1 - share);
            hist[idx1] += _increment * share;
        }
        else
        {
            hist[idx1] += _increment;
        }
    }
    return hist;
}
    
//the first vector corresponds to the first image
//to efficiently compute the histogram for the second image
vector<double> MutualInformation::computeHist2d(const vector<double> & valVec1, const vector<double> & valVec2) const
{
    //row-major 2D grid
    assert(valVec1.size() == valVec2.size());
    vector<double> hist(_numBins * _numBins, 0);
    for (int i = 0; i < valVec1.size(); i++)
    {
        int idx11, idx21;
        int idx12, idx22;
        double share1, share2;
        computeShares(valVec1[i], idx11, idx12, share1);
        computeShares(valVec2[i], idx21, idx22, share2);
        if (idx12 != -1 and idx22 != -1) 
        {
            hist[idx21 * _numBins + idx11] += _increment * share1*share2;
            hist[idx21 * _numBins + idx12] += _increment *(1 - share1)*share2;
            hist[idx22 * _numBins + idx11] += _increment *(1 - share2)*share1;
            hist[idx22 * _numBins + idx12] += _increment *(1 - share1)*(1 - share2);
        }
        else if (idx12 != -1) 
        {
            hist[idx21 * _numBins + idx11] += share1*_increment;
            hist[idx21 * _numBins + idx12] += (1 - share1)*_increment;
        }
        else if (idx22 != -1) 
        {
            hist[idx21 * _numBins + idx11] += _increment*share2;
            hist[idx22 * _numBins + idx11] += (1 - share2)*_increment;
        }
        else
        {
            hist[idx21 * _numBins + idx11] += _increment;
        }
    }
    return hist;
}

vector<double> MutualInformation::reduceHist(const vector<double> & hist2d) const
{
    assert(hist2d.size() % _numBins == 0);
    vector<double> hist;
    hist.reserve(_numBins);
    int i = 0;
    for (auto iter = hist2d.begin(); iter != hist2d.end(); iter += _numBins)
    {
        hist.push_back(accumulate(iter, iter +  _numBins, 0.));
    }
    return hist;
}

