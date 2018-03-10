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

#include "localization/cost_function_mi.h"

#include "std.h"
#include "eigen.h"
#include "ocv.h" //TODO replace here Mat32f by a pointer
#include "ceres.h"
#include "io.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "projection/jacobian.h"
#include "reconstruction/triangulator.h"

bool MutualInformation::Evaluate(double const * parameters,
        double * cost, double * gradient) const
{
    const int POINT_NUMBER = _dataPack.cloud.size();
    for (int i = 0; i < 6; i++)
    {
        if (std::isnan(parameters[i]) or std::isinf(parameters[i])) return false;
    }
    Transf xiBase(parameters);
    Transf xiCam = xiBase.compose(_xiBaseCam);
    
    // point cloud in frame 2
    vector<Vector3d> transformedPoints;
    xiCam.inverseTransform(_dataPack.cloud, transformedPoints);
    // init the image interpolation
    ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(_imageGrid);
    
    bool computeGrad = (gradient != NULL);
    
    vector<double> valVec2(POINT_NUMBER, 0);
    vector<Covector2d> gradVec;
    if (computeGrad)
    {
        gradVec.resize(POINT_NUMBER, Covector2d(0, 0));
    }
    *cost = 0;
    for (int i = 0; i < POINT_NUMBER; i++)
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
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale,
                        &f, &grad[1], &grad[0]);
                grad *= _invScale;  // normalize according to the scale
            }
            else 
            {
                double & f = valVec2[i];
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale, &f);
            }
            
        }
        
    }
    vector<double> hist12 = computeHist2d(_dataPack.valVec, valVec2);
    vector<double> hist2 = reduceHist(hist12);
    vector<double> logVec12(_numBins * _numBins, 0);
    
    // compute the cost
    *cost = 0;
    for (int idx2 = 0; idx2 < _numBins; idx2++)
    {
        for (int idx1 = 0; idx1 < _numBins; idx1++)
        {
            const double & p12 = hist12[idx2 * _numBins + idx1];
            const double log12 = log(p12 / (hist2[idx2] * _hist1[idx1]));
            if (p12 > 0)
            {
                logVec12[idx2 * _numBins + idx1] = log12;
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
        CameraJacobian jacobianCalculator(_camera, xiBase, _xiBaseCam);
        for (int i = 0; i < POINT_NUMBER; i++)
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
                    dMIdP = logVec12[idx21 * _numBins + idx11] * share1
                        + logVec12[idx21 * _numBins + idx12] * (1 - share1)
                        - logVec12[idx22 * _numBins + idx11] * share1
                        - logVec12[idx22 * _numBins + idx12] * (1 - share1);
                }
                else
                {
                    dMIdP = logVec12[idx21 * _numBins + idx11]
                        - logVec12[idx22 * _numBins + idx11];
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


MutualInformationOdom::MutualInformationOdom(  const ICamera * camera, 
                        const PhotometricPack & dataPack, 
                        const Transf xiBaseCam,
                        const Transf xiOdom,
                        const Transf xiPrior,
                        const Mat32f & img2, 
                        double scale, 
                        int numBins, 
                        double valMax,
                        const double errV, 
                        const double errW,
                        const double lambdaT, 
                        const double lambdaR
            ) :
            MutualInformation(camera, dataPack, xiBaseCam, img2, scale, numBins, valMax),
            _xiPrior(xiPrior),
            _C(Matrix6d::Zero())
{
    const double delta = xiOdom.rot()(2);
    const double l = xiOdom.trans().norm();
    
    const double delta2 = delta / 2.;
    const double l2 = l / 2.;
    
    const double s = sin(delta2);
    const double c = cos(delta2);
    
    Matrixd<3, 2> dfdu;

    dfdu <<     s,      -l2 * c,
                c,      l2 * s,                   
                0,           1;

    Matrix2d Cu;
    Cu <<   errV * errV * l * l,              0,
                      0,    errW * errW * delta * delta; 
   
    Matrix3d lambdaMat = Matrix3d::Identity();
    lambdaMat(0, 0) *= lambdaT * lambdaT;
    lambdaMat(1, 1) *= lambdaT * lambdaT;
    lambdaMat(2, 2) *= lambdaR * lambdaR;
    Matrix3d Cx = dfdu * Cu * dfdu.transpose() + lambdaMat;
    Matrix3d CxInv = Cx.inverse();
    _C.topLeftCorner<2, 2>() = CxInv.topLeftCorner<2, 2>();
    _C.topRightCorner<2, 1>() = CxInv.topRightCorner<2, 1>();
    _C.bottomLeftCorner<1, 2>() = CxInv.bottomLeftCorner<1, 2>();
    _C(5, 5) = CxInv(2, 2);
    _C(2, 2) = 1 / (lambdaT * lambdaT);
    _C(3, 3) = 1 / (lambdaR * lambdaR);
    _C(4, 4) = 1 / (lambdaR * lambdaR);
    
    Matrix3d M = interOmegaRot(xiPrior.rot());
    Matrix3d R = xiPrior.rotMatInv();
    
    _J.topLeftCorner<3, 3>() = _C.topLeftCorner<3, 3>() * R;
    _J.topRightCorner<3, 3>() = _C.topRightCorner<3, 3>() * R * M;
    _J.bottomLeftCorner<3, 3>() = _C.bottomLeftCorner<3, 3>() * R;
    _J.bottomRightCorner<3, 3>() = _C.bottomRightCorner<3, 3>() * R * M;
    _C *= 0.5; //to avoid division by 2 at every evaluation
}

bool MutualInformationOdom::Evaluate(double const * parameters, double * cost, double * gradient) const
{
    MutualInformation::Evaluate(parameters, cost, gradient);
    Transf xi(parameters);
    Covector6d err;
    _xiPrior.inverseCompose(xi).toArray(err.data());
    
    const double DAMPING = 0.0002;
    
    Vector6d priorGrad = err * _J;
    for (int i = 0; i < 6; i++)
    {
        gradient[i] += priorGrad[i] * DAMPING;
    }
    *cost += DAMPING * double(err * _C * err.transpose());
}


