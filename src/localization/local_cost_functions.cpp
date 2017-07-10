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

#include "localization/local_cost_functions.h"

#include "std.h"
#include "eigen.h"
#include "ocv.h" //TODO replace here Mat32f by a pointer
#include "ceres.h"
#include "io.h"

#include "geometry/geometry.h"
#include "projection/generic_camera.h"
#include "projection/jacobian.h"
#include "reconstruction/triangulator.h"

PhotometricCostFunction::PhotometricCostFunction(const ICamera * camera, const Transf & xiBaseCam,
            const PhotometricPack & dataPack,
            const Mat32f & img2, double scale) :
            _camera(camera->clone()),
            _dataPack(dataPack),
            _xiBaseCam(xiBaseCam),
            _imageGrid(img2.cols, img2.rows, (float*)(img2.data)),
            _invScale(1. / scale) 
    {
        mutable_parameter_block_sizes()->clear();
        mutable_parameter_block_sizes()->push_back(6);
        set_num_residuals(_dataPack.cloud.size());
    }
    
void PhotometricCostFunction::lossFunction(const double x, double & rho, double & drhodx) const
{
    double s = sign(x);
    const double e = exp(-abs(x) / LOSS_FACTOR);
    rho = s * LOSS_FACTOR * (1. - e);
    drhodx = e;

    return;
}

/*
A cost function with analytic jacobian
works faster than autodiff version and works with any ICamera
*/
    
bool PhotometricCostFunction::Evaluate(double const * const * parameters,
        double * residual, double ** jacobian) const
{
    const int POINT_NUMBER = _dataPack.cloud.size();
    const double NORMALIZER = 1.;
    Transf xiBase(parameters[0]);
    Transf xiCam = xiBase.compose(_xiBaseCam);
    // point cloud in frame 2
    vector<Vector3d> transformedPoints;
    xiCam.inverseTransform(_dataPack.cloud, transformedPoints);
    
    // init the image interpolation
    ceres::BiCubicInterpolator<Grid2D<float>> imageInterpolator(_imageGrid);
    
    bool computeJac = (jacobian != NULL and jacobian[0] != NULL);
    if (computeJac)
    {
        // L_uTheta
        CameraJacobian jacobianCalculator(_camera, xiBase, _xiBaseCam);
        for (int i = 0; i < POINT_NUMBER; i++)
        {
            Vector2d pt;
            if (_camera->projectPoint(transformedPoints[i], pt)) 
            {
                double f;
                // image interpolation and gradient
                Covector2d grad;
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale,
                        &f, &grad[1], &grad[0]);
                grad *= _invScale;  // normalize according to the scale
                
                residual[i] = (f - _dataPack.valVec[i]); // / (_dataPack.valVec[i] + 10);
//                if (abs(residual[i]) > 10)
//                {
//                    residual[i] = 0.;
//                    fill(jacobian[0] + i*6, jacobian[0] + i*6 + 6, 0.);
//                }
//                else
//                {                
                jacobianCalculator.dfdxi(transformedPoints[i], grad, jacobian[0] + i*6);
                double k;
                lossFunction(residual[i], residual[i], k);
                residual[i] *= NORMALIZER;
                k *= NORMALIZER;
                for (int j = 0; j < 6; j++)
                {
                    jacobian[0][i*6 + j] *= k;
                }
//                }
                
//                for (int j = 0; j < 6; j++) *(jacobian[0] + i*6 + j) /= (_dataPack.valVec[i] + 10); 
                
                
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
        for (int i = 0; i < POINT_NUMBER; i++)
        {
            Vector2d pt;
            if (_camera->projectPoint(transformedPoints[i], pt)) 
            {
                double f;
                imageInterpolator.Evaluate(pt[1] * _invScale, pt[0] * _invScale, &f);
                residual[i] = (f - _dataPack.valVec[i]);// / (_dataPack.valVec[i] + 10);
                double k;
                lossFunction(residual[i], residual[i], k);
                residual[i] *= NORMALIZER;
//                if (abs(residual[i]) > 10) residual[i] = 0.;
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

bool MonoReprojectCost::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //compute the transformation chain
    Transf xiOdom(params[0]);
    Transf xi21 = _xiBaseCam.inverseCompose(xiOdom.inverseCompose(_xiBaseCam));
    
    //compute the points in the camera frame
    Vector3dVec xVec2 =_xVec1;
    for (int i = 0; i < 5; i++)
    {
        xVec2[i] *= params[1][i];
    }
    xi21.transform(xVec2, xVec2);
    
    //compute the reprojection error
    for (int i = 0; i < 5; i++)
    {
        Vector2d modProj;
        if (_camera->projectPoint(xVec2[i], modProj)) 
        {
            Map<Vector2d> diff(residual + 2*i);
            diff = modProj - _pVec2[i];
        }
        else
        {
            residual[2*i + 1] = residual[2*i] = DOUBLE_BIG;
        }
    }
    
    if (jacobian != NULL)
    {
        //odometry jacobian
        if (jacobian[0] != NULL) 
        {
            InterJacobian jacobianCalculator(_camera, _xiBaseCam.inverse(),
                                    xiOdom, JAC_INVERTED);
            
            for (int i = 0; i < 5; i++)
            {
                jacobianCalculator.dpdxi(xVec2[i], jacobian[0] + i*12,
                                                   jacobian[0] + i*12 + 6);
            }
        }
        
        //length jacobian
        if (jacobian[1] != NULL)
        {
            Matrix3d R21 = xi21.rotMat();
            fill(jacobian[1], jacobian[1] + 50, 0);
            for (int i = 0; i < 5; i++)
            {
                Matrix23drm dpdx;
                Vector3d n2 = R21 * _xVec1[i];
                _camera->projectionJacobian(xVec2[i], dpdx.data(), dpdx.data() + 3);
                Vector2d dpdl = dpdx * n2;
                jacobian[1][i*11] = dpdl[0];
                jacobian[1][i*11 + 5] = dpdl[1];
            }
        }
    }
    return true;
}


bool SparseReprojectCost::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //compute the transformation chain
    Transf xiOdom(params[0]);
    Transf xi12 = _xiBaseCam.inverseCompose(xiOdom.compose(_xiBaseCam));
    
    //compute the points in the camera frame
    vector<double> lambdaVec(_xVec1.size());
    vector<double> jacVec;
    Triangulator triangulator(xi12);
    if (jacobian != NULL and jacobian[0] != NULL) 
    {
        jacVec.resize(_xVec1.size() * 6);
        triangulator.computeRegular(_xVec1, _xVec2, lambdaVec.data(), NULL, jacVec.data(), NULL);
    }
    else
    {
        triangulator.computeRegular(_xVec1, _xVec2, lambdaVec.data());
    }
    
    Vector3dVec xVec2 =_xVec1;
    for (int i = 0; i < xVec2.size(); i++)
    {
        xVec2[i] *= lambdaVec[i];
    }
    xi12.inverseTransform(xVec2, xVec2);
    
    //compute the reprojection error
    for (int i = 0; i < xVec2.size(); i++)
    {
        Vector2d modProj;
        if (_camera->projectPoint(xVec2[i], modProj)) 
        {
            Map<Vector2d> diff(residual + 2*i);
            diff = modProj - _pVec2[i];
        }
        else
        {
            residual[2*i + 1] = residual[2*i] = DOUBLE_BIG;
        }
    }
    
    if (jacobian != NULL)
    {
        
        if (jacobian[0] != NULL) 
        {
            //odometry jacobian
            InterJacobian jacobianCalculator(_camera, _xiBaseCam.inverse(),
                                    xiOdom, JAC_INVERTED);
            
            for (int i = 0; i < xVec2.size(); i++)
            {
                if (residual[2*i] == DOUBLE_BIG)
                {
                    fill(jacobian[0] + i * 12, jacobian[0] + (i + 1) * 12, 0);
                }
                else
                {
                    jacobianCalculator.dpdxi(xVec2[i], jacobian[0] + i*12,
                                                   jacobian[0] + i*12 + 6);
                }
            }
            
            //length jacobian
            Matrix3d R21 = xi12.rotMatInv();
            Matrix3d RcamBase = _xiBaseCam.rotMatInv();
            
            //jacobian given by triangulate is computed wrt ( v_1_2 | om_1_2 )
            
            // om_1_2 = M * r_dot
            Matrix3d M = RcamBase * interOmegaRot(xiOdom.rot());    
            
            // v_2 = v_b + om_0_b x t_b_c
            // v_1_2 = R_c_b * t_dot + Q * r_dot
            Vector3d tBaseCam1 = RcamBase * R21.transpose() * _xiBaseCam.trans();
            Matrix3d Q = -hat(tBaseCam1) * M;
            for (int i = 0; i < xVec2.size(); i++)
            {
                if (residual[2*i] == DOUBLE_BIG) continue;
                Matrix23drm dpdx;
                Vector3d n2 = R21 * _xVec1[i];
                _camera->projectionJacobian(xVec2[i], dpdx.data(), dpdx.data() + 3);
                Vector2d dpdl = dpdx * n2;
                Map<Covector3d> dudt(jacobian[0] + i*12);
                Map<Covector3d> dudr(jacobian[0] + i*12 + 3);
                Map<Covector3d> dvdt(jacobian[0] + i*12 + 6);
                Map<Covector3d> dvdr(jacobian[0] + i*12 + 9);
                
                Map<Covector3d> dldv(jacVec.data() + i*6);
                Map<Covector3d> dldw(jacVec.data() + i*6 + 3);
                Covector3d dldt = dldv * RcamBase;
                Covector3d dldr = dldw * M + dldv * Q;
                dudt += dpdl[0] * dldt;
                dudr += dpdl[0] * dldr;
                dvdt += dpdl[1] * dldt;
                dvdr += dpdl[1] * dldr;
            }
        }
    }
    return true;
}

OdometryPrior::OdometryPrior(const double errV, const double errW,
        const double lambdaT, const double lambdaR,
        const Transf xiOdom) : 
    _A(Matrix6d::Zero()),
    _xiPrior(xiOdom)
{
    const double delta = xiOdom.rot().norm();
    const double l = xiOdom.trans().norm();
    
    const double delta2 = delta / 2.;
    const double l2 = l / 2.;
    
    const double s = sin(delta2);
    const double c = cos(delta2);
    
    Matrixd<3, 2> dfdu;

    dfdu <<     c,     l2 * s,  
                -s,      l2 * c, 
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
    Eigen::LLT<Matrix3d> lltOfCxInv(CxInv); // compute the Cholesky decomposition of A
    Matrix3d U = lltOfCxInv.matrixU();
//    cout << U.transpose() * U << endl;
//    cout << CxInv << endl;

    _A.topLeftCorner<2, 2>() = U.topLeftCorner<2, 2>();
    _A.topRightCorner<2, 1>() = U.topRightCorner<2, 1>();
    _A(2, 2) = 1. / lambdaT;
    _A(3, 3) = 1. / lambdaR;
    _A(4, 4) = 1. / lambdaR;
    _A(5, 5) = U(2, 2);
    
    Matrix3d M = interOmegaRot(xiOdom.rot());
    Matrix3d R = xiOdom.rotMatInv();
    
    _J.topLeftCorner<3, 3>() = _A.topLeftCorner<3, 3>() * R;
    _J.topRightCorner<3, 3>() = _A.topRightCorner<3, 3>() * R * M;
    _J.bottomLeftCorner<3, 3>() = Matrix3d::Zero();
    _J.bottomRightCorner<3, 3>() = _A.bottomRightCorner<3, 3>() * R * M;
//    cout << U << endl;
//    cout << _A << endl;
//        assert(false);
}

bool OdometryPrior::Evaluate(double const * const * params,
        double * residual, double ** jacobian) const
{
    //TODO replace the difference between two transformation by InverseCompose
    Transf xi(*params);
    Transf delta = _xiPrior.inverseCompose(xi);
    Vector6d vecDelta;
    delta.toArray(vecDelta.data());
    Map<Vector6d> res(residual);
    res = _A * vecDelta;
    
    if (jacobian != NULL)
    {
        //first transformation
        if (jacobian[0] != NULL)
        {
            Map<Matrix6drm> jac(jacobian[0]);
            jac = _J;
        }
    }
    return true;    
}
