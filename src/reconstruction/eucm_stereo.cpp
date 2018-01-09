/*
This file is part of visgeom.

visgeom is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

visgeom is distributed in the hope that it will be useful,
but WITHOUdouble ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with visgeom.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#include "reconstruction/eucm_stereo.h"

StereoParameters::StereoParameters(const ptree & params) :
        ScaleParameters(params)
{
    for (auto & item : params.get_child("stereo_parameters"))
    {
        const string & pname = item.first;
        if (pname == "disparity_max")       dispMax = item.second.get_value<int>();
        else if (pname == "error_max")      maxError = item.second.get_value<int>();
        else if (pname == "verbosity")      verbosity = item.second.get_value<int>();
        else if (pname == "hypotheses")     hypMax = item.second.get_value<int>();
        else if (pname == "hypo_difference")    maxHypDiff = item.second.get_value<int>();
        else if (pname == "flaw_cost")    flawCost = item.second.get_value<int>();
        else if (pname == "descriptor_size")    descLength = item.second.get_value<int>();
        else if (pname == "scales")    scaleVec = readVector<int>(item.second);
        else if (pname == "descriptor_response_thresh") descRespThresh = item.second.get_value<int>();
        else if (pname == "num_epipolar_planes") numEpipolarPlanes = item.second.get_value<int>();
        else if (pname == "epipole_margin") epipoleMargin = pow(item.second.get_value<int>(), 2);
    }
}

EnhancedStereo::EnhancedStereo(const EnhancedCamera * cam1, const EnhancedCamera * cam2,
        const StereoParameters & params) :
    _params(params),
    _camera1(cam1->clone()),
    _camera2(cam2->clone()),
    _epipolarCurves(cam1, cam2, _params.numEpipolarPlanes, params.verbosity),
    _epipolarDescriptor(params.descLength, params.descRespThresh, params.scaleVec),
    HALF_LENGTH(params.descLength / 2),
    MARGIN(params.descLength - 1),
    _triangulator(1e-3)
{ 
    assert(params.descLength % 2 == 1);
}

EnhancedStereo::~EnhancedStereo()
{
    delete _camera1;
    _camera1 = NULL;
    delete _camera2;
    _camera2 = NULL;
}


void EnhancedStereo::setTransformation(const Transf & T12) 
{
    _transf12 = T12; 
    _R12 = T12.rotMat();
    _R21 = T12.rotMatInv();
    _t12 = T12.trans();
    _triangulator.setTransformation(T12);
    _epipolarCurves.setTransformation(T12);
}

int computeError(int v, int thMin, int thMax)
{
    return max(0, max(thMin - v, v - thMax));
}

vector<int> compareDescriptor(const vector<uint8_t> & desc,
        const vector<uint8_t> & sampleVec, int flawCost)
{
    vector<int> thMinVec(desc.size()), thMaxVec(desc.size());
    
    for (int i = 1; i < desc.size() - 1; i++)
    {
        const int d = desc[i];
        int d1 = (desc[i] + desc[i - 1]) / 2;
        int d2 = (desc[i] + desc[i + 1]) / 2;
        thMinVec[i] = min(d, min(d1, d2));
        thMaxVec[i] = max(d, max(d1, d2));
    }
    
    if (desc[0] > desc[1])
    {
        thMinVec[0] = (desc[0] + desc[1]) / 2;
        thMaxVec[0] = desc[0]; 
    }
    else
    {
        thMaxVec[0] = (desc[0] + desc[1]) / 2;
        thMinVec[0] = desc[0];
    }
    
    const int d = desc[desc.size() - 2];
    if (desc.back() > d)
    {
        thMinVec.back() = (desc.back() + d) / 2;
        thMaxVec.back() = desc.back(); 
    }
    else
    {
        thMaxVec.back() = (desc.back() + d) / 2;
        thMinVec.back() = desc.back();
    }
    
    
    const int HALF_LENGTH = desc.size() / 2;
    vector<int> rowA(sampleVec.size()), rowB(sampleVec.size());
    /////////////////////////////////////////////////////////////////
    /*
    //match the first half
    for (int i = 0; i < sampleVec.size(); i++)
    {
        rowA[i] = abs(int(sampleVec[i]) - int(desc[0]));
    }
    for (int i = 1; i <= HALF_LENGTH; i++)
    {
        rowB[0] = rowA[0] + flawCost + abs(int(sampleVec[0]) - int(desc[i]));
        int cost = min(rowA[1] + flawCost, rowA[0]);
        rowB[1] = cost + abs(int(sampleVec[1]) - int(desc[i]));
        for (int j = 2; j < sampleVec.size(); j++)
        {
            cost = min(min(rowA[j] + flawCost, rowA[j - 1]), rowA[j - 2] + flawCost);
            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        }
        swap(rowA, rowB);
    }
    vector<int> rowC(sampleVec.size()); //center cost
    swap(rowA, rowC);
    
    //match the second half (from the last pixel to first)
    for (int i = 0; i < sampleVec.size(); i++)
    {
        rowA[i] = abs(sampleVec[i] - desc.back());
    }
    for (int i = desc.size() - 2; i > HALF_LENGTH; i--)
    {
        for (int j = 0; j < sampleVec.size() - 2; j++)
        {
            int cost = min(min(rowA[j] + flawCost, rowA[j + 1]), rowA[j + 2] + flawCost);
            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        }
        int j = sampleVec.size() - 2;
        int cost = min(rowA[j] + flawCost, rowA[j + 1]);
        rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        rowB.back() = rowA.back() + flawCost + abs(int(sampleVec.back()) - int(desc[i]));
        swap(rowA, rowB);
    }
    */
    ////////////////////////////////////////////////////////////////////////////////
    
        //match the first half
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(int(sampleVec[i]) - int(desc[0]));
        rowA[i] = computeError(sampleVec[i], thMinVec[0], thMaxVec[0]);
    }
    for (int i = 1; i <= HALF_LENGTH; i++)
    {
//        rowB[0] = rowA[0] + flawCost + abs(int(sampleVec[0]) - int(desc[i]));
        rowB[0] = rowA[0] + flawCost + computeError(sampleVec[0], thMinVec[i], thMaxVec[i]);
        int cost = min(rowA[1] + flawCost, rowA[0]);
//        rowB[1] = cost + abs(int(sampleVec[1]) - int(desc[i]));
        rowB[1] = cost + computeError(sampleVec[1], thMinVec[i], thMaxVec[i]);
        for (int j = 2; j < sampleVec.size(); j++)
        {
            cost = min(min(rowA[j] + flawCost, rowA[j - 1]), rowA[j - 2] + flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        swap(rowA, rowB);
    }
    vector<int> rowC(sampleVec.size()); //center cost
    swap(rowA, rowC);
    
    //match the second half (from the last pixel to first)
    for (int i = 0; i < sampleVec.size(); i++)
    {
//        rowA[i] = abs(sampleVec[i] - desc.back());
        rowA[i] = computeError(sampleVec[i], thMinVec.back(), thMaxVec.back());
    }
    for (int i = desc.size() - 2; i > HALF_LENGTH; i--)
    {
        for (int j = 0; j < sampleVec.size() - 2; j++)
        {
            int cost = min(min(rowA[j] + flawCost, rowA[j + 1]), rowA[j + 2] + flawCost);
//            rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
            rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        }
        int j = sampleVec.size() - 2;
        int cost = min(rowA[j] + flawCost, rowA[j + 1]);
//        rowB[j] = cost + abs(int(sampleVec[j]) - int(desc[i]));
        rowB[j] = cost + computeError(sampleVec[j], thMinVec[i], thMaxVec[i]);
        
//        rowB.back() = rowA.back() + flawCost + abs(int(sampleVec.back()) - int(desc[i]));
        rowB.back() = rowA.back() + flawCost + computeError(sampleVec.back(), thMinVec[i], thMaxVec[i]);
        swap(rowA, rowB);
    }
    
    //accumulate the cost
    for (int i = 0; i < sampleVec.size() - 2; i++)
    {
        rowC[i] += min(min(rowA[i] + flawCost, rowA[i + 1]), rowA[i + 2] + flawCost);
    }
    int i = rowC.size() - 2;
    rowC[i] += min(rowA[i] + flawCost, rowA[i + 1]);
    rowC.back() += rowA.back() + flawCost;
    return rowC;
} 


double EnhancedStereo::triangulate(double u1, double v1, double u2,
        double v2, CameraIdx camIdx) const
{
    if (_params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    Vector3d p, q;
    if (not _camera1->reconstructPoint(Vector2d(u1, v1), p) or 
        not _camera2->reconstructPoint(Vector2d(u2, v2), q) )
    {
        if (_params.verbosity > 2) 
        {
            cout << "    not reconstructed " << Vector2d(u1, v1).transpose(); 
            cout << " # " << Vector2d(u2, v2).transpose() << endl;
        }
        return OUT_OF_RANGE;
    }
    double lambda;
    if (camIdx == CAMERA_1)
    {
        _triangulator.computeRegular(p, q, &lambda);
        return p.norm() * lambda; //TODO change to lambda only?
    }
    else
    {
        _triangulator.computeRegular(p, q, NULL, &lambda);
        return q.norm() * lambda; //TODO change to lambda only?
    }
}

//FIXME
const double TRIANGULATE_DIST_MAX = 100;
bool EnhancedStereo::triangulate(const double u1, const double v1, const double u21, const double v21,
        const double u22, const double v22, double & d, double & sigma, CameraIdx camIdx) const
{
    if (_params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    Vector3d p, q1, q2;
    if (not _camera1->reconstructPoint(Vector2d(u1, v1), p) or 
        not _camera2->reconstructPoint(Vector2d(u21, v21), q1) or 
        not _camera2->reconstructPoint(Vector2d(u22, v22), q2))
    {
        if (_params.verbosity > 2) 
        {
            cout << "    not reconstructed " << Vector2d(u1, v1).transpose(); 
            cout << " # " << Vector2d(u21, v21).transpose() << endl;
        }
        return false;
    }
    
    double lambda1, lambda2;
    if (camIdx == CAMERA_1)
    {
        _triangulator.computeRegular(p, q1, &lambda1);
        _triangulator.computeRegular(p, q2, &lambda2);
        double pnorm = p.norm();
        lambda1 *= pnorm; //TODO change to lambda only?
        lambda2 *= pnorm; 
    }
    else
    {
        _triangulator.computeRegular(p, q1, NULL, &lambda1);
        _triangulator.computeRegular(p, q2, NULL, &lambda2);
        lambda1 *= q1.norm(); //TODO change to lambda only?
        lambda2 *= q2.norm(); 
    }
    
    //result
    if (lambda1 < TRIANGULATE_DIST_MAX) 
    {
        sigma = abs(lambda2 - lambda1);
        d = lambda1;
    }
    else
    {
        sigma = OUT_OF_RANGE;
        d = OUT_OF_RANGE;
    }
}

bool EnhancedStereo::triangulate(const Vector2d pt11, const Vector2d pt12, const Vector2d pt21,
        const Vector2d pt22, double & d, double & sigma) const
{
    if (_params.verbosity > 3) cout << "EnhancedStereo::triangulate" << endl;
    Vector3d p1, p2, q1, q2;  //spatial direction points
    if (not _camera1->reconstructPoint(pt11, p1) or 
        not _camera1->reconstructPoint(pt12, p2) or
        not _camera2->reconstructPoint(pt21, q1) or 
        not _camera2->reconstructPoint(pt22, q2))
    {
        if (_params.verbosity > 2) 
        {
            cout << "    not reconstructed " << pt11.transpose(); 
            cout << " # " << pt21.transpose() << endl;
        }
        return false;
    }
    
    double lambda1, lambda2;
    _triangulator.computeRegular(p1, q1, &lambda1);
    _triangulator.computeRegular(p2, q2, &lambda2);
    double pnorm = p1.norm();
    lambda1 *= pnorm; //TODO change to lambda only?
    lambda2 *= pnorm; 
    
    //result
    sigma = abs(lambda2 - lambda1);
    d = lambda1;
}
