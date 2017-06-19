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
Depth container. 
NOTE:
(u, v) is an image point 
(x, y) is a depth map point
*/

#include "reconstruction/depth_map.h"

#include "io.h"
#include "std.h"
#include "eigen.h"

//FIXME remove hypStep

void DepthMap::applyMask(const Mat8u & mask)
{
    for (int y = 0; y < yMax; y++)
    {
        for (int x = 0; x < xMax; x++)
        {
            if ( not mask(vConv(y), uConv(x)) ) 
            {
                for (int h = 0; h < hMax; h++)
                {
                    at(x, y, h) = 0;
                }
            }
        }
    }
}

//check the limits
bool DepthMap::isValid(const int x, const int y, const int h) const
{
    return ( (x >= 0) and (x < xMax) and (y >= 0) and (y < yMax) and (h >= 0) and (h < hMax) );
}

bool DepthMap::isValid(const Vector2d pt, const int h) const
{
    return isValid(pt[0], pt[1], h);
}

void DepthMap::swapHypotheses(const int x, const int y, const int h1, const int h2)
{
    const double temp_d = at(x, y, h1);
    const double temp_sigma = sigma(x, y, h1);
    const double temp_cost = cost(x, y, h1);
    at(x, y, h1) = at(x, y, h2);
    sigma(x, y, h1) = sigma(x, y, h2);
    cost(x, y, h1) = cost(x, y, h2);
    at(x, y, h2) = temp_d;
    sigma(x, y, h2) = temp_sigma;
    cost(x, y, h2) = temp_cost;
}

void DepthMap::sortHypStack(const int x, const int y)
{
    for (int h1 = 0; h1 < hMax - 1; ++h1)
    {
        for (int h2 = h1 + 1; h2 < hMax; ++h2)
        {
            //Swap if previous element is default and next element is valid hyp
            if ( (at(x, y, h1) <= MIN_DEPTH) && (at(x, y, h2) > MIN_DEPTH) )
            {
                swapHypotheses(x, y, h1, h2);
            }
            //Swap to order in ascending order of cost
            else if ( cost(x, y, h1) > cost(x, y, h2) )
            {
                swapHypotheses(x, y, h1, h2);
            }
        }
    }
}

bool DepthMap::pushHypothesis(const Vector3d & X, const double sigmaVal)
{
    Vector2d pt;
    if (not cameraPtr->projectPoint(X, pt)) return false;
    
    const int x = xConv(pt[0]);
    const int y = yConv(pt[1]);
    const double d = X.norm();

    return pushHypothesis(x, y, d, sigmaVal);
}

bool DepthMap::pushHypothesis(const int x, const int y, const double d, const double sigmaVal)
{
    if (not isValid(x, y)) return false;
    int h = 0;
    // while (h < hMax)
    // {
    //     if ( (at(x, y, h) > MIN_DEPTH) and (cost(x, y, h) <= DEFAULT_COST_DEPTH) ) h++;
    //     else break;
    // }
    // if (h == hMax) return false;    

    //Insert element into stack, so that stack is always sorted in cost ascending order
    if ( (at(x, y, hMax-1) > MIN_DEPTH) and (cost(x, y, hMax-1) <= DEFAULT_COST_DEPTH) ) return false;
    h = hMax-1; //Replace element at the end of the stack
    at(x, y, h) = d;
    sigma(x, y, h) = sigmaVal;
    cost(x, y, h) = DEFAULT_COST_DEPTH + 2*COST_CHANGE;
    sortHypStack(x, y);
    // double temp_d = d;
    // double temp_sigma = sigmaVal;
    // double temp_cost = DEFAULT_COST_DEPTH;
    // double temp2_d;
    // double temp2_sigma;
    // double temp2_cost;
    // while (h < hMax)
    // {
    //     temp2_d = at(x, y, h);
    //     temp2_sigma = sigma(x, y, h);
    //     temp2_cost = cost(x, y, h);
    //     at(x, y, h) = temp_d;
    //     sigma(x, y, h) = temp_sigma;
    //     cost(x, y, h) = temp_cost;
    //     temp_d = temp2_d;
    //     temp_sigma = temp2_sigma;
    //     temp_cost = temp2_cost;
    //     h++;
    // }
    return true;
}

bool DepthMap::pushImageHypothesis(const int u, const int v, const double d, const double sigmaVal)
{
    return pushHypothesis(xConv(u), yConv(v), d, sigmaVal);
}

bool DepthMap::match(const double v1, const double s1, const double v2, const double s2)
{
    const double delta = abs(v1 - v2);
    const bool output = delta < 1 * s1 or delta < 1 * s2;
    // if (output and ((v1<1) or (v2<1)))
    // {
    //     cout << "v1: " << v1 << " v2: " << v2 << " s1: " << s1 << " s2: " << s2 << endl;
    // }
    return output;
}

void DepthMap::filter(double & v1, double & s1, const double v2, const double s2)
{
    double denom = s1 + s2;
    v1 = (v1 * s2 + v2 * s1) / denom;
    s1 = s1 * s2 / denom;
}

bool DepthMap::filterPushHypothesis(const Vector3d & X, const double sigmaVal)
{
    Vector2d pt;
    if (not cameraPtr->projectPoint(X,pt)) return false;

    const int x = xConv(pt[0]);
    const int y = yConv(pt[1]);
    const double d = X.norm();

    return filterPushHypothesis(x, y, d, sigmaVal);
}

bool DepthMap::filterPushHypothesis(const int x, const int y, const double d, const double sigmaVal)
{
    if (not isValid(x, y)) return false;
    bool matchFound = false;
    bool hypExist = false;
    for(int h = 0; h < hMax; ++h)
    {
        double & d1 = at(x, y, h);
        double & sigma1 = sigma(x, y, h);
        double & cost1 = cost(x, y, h);
        if( d1 == OUT_OF_RANGE )
        {
            if ( h == 0 )
            {
                pushHypothesis(x, y, d, sigmaVal);
    // if(d<1) cout << "Pushing: x: " << x << " y: " << y << " d: " << d << " s: " << sigmaVal << endl;
                return true;
            }
            break;
        }
        else if( (not matchFound) and match(d1, sigma1, d, sigmaVal) )
        {
            filter(d1, sigma1, d, sigmaVal);
            // if(d1<1)
            // {
            //     cout << "After merge: d1: " << d1 << "  s1: " << sigma1 << endl;
            // }
            // cost1 = max(cost1 - 1.0, 0.0);
            cost1 -= 2*COST_CHANGE;
            cost1 = std::max(cost1, 0.0);
            matchFound = true;
            // std::cout << "Merged " << x << " " << y << " " << h << " " << cost1 << endl;
        }
        else
        {
            cost1 += COST_CHANGE; // increase cost for all hypotheses that did not match
        }
    }
    sortHypStack(x, y);
    // if (matchFound) return true;
    //The program reaches this area if no matches were found
    // else return pushHypothesis(x, y, d, sigmaVal);
    return matchFound;
    // else 
    // {
    //     // cout << "Mismatch: " << x << " " << y << " " << d << " " << sigmaVal << endl;
    //     // cout << "\t Options:" << endl;
    //     // for(int h = 0; h < hMax; ++h)
    //     //     cout << "\t" << at(x, y, h) << " " << sigma(x, y, h) << " " << cost(x, y, h) << endl;
    //     return false;
    // }
}

void DepthMap::costRejection(const double rejectionThreshold)
{
    // int count = 0;
    for (int h = 0; h < hMax; ++h)
    {
        for (int y = 0; y < yMax; ++y)
        {
            for (int x = 0; x < xMax; ++x)
            {
                double & c = cost(x, y, h);
                // c += costChange;
                if( c > rejectionThreshold )
                {
                    at(x, y, h) = DEFAULT_DEPTH;
                    sigma(x, y, h) = DEFAULT_SIGMA_DEPTH;
                    c = DEFAULT_COST_DEPTH;
                    // count++;
                }
            }
        }
    }
    // cout << "discarded " << count << " hyps" << endl;
}

// nearest neighbor interpolation
double DepthMap::nearest(const int u, const int v, const int h) const
{
    int xd = xConv(u);
    int yd = yConv(v);
    if (isValid(xd, yd, h)) return at(xd, yd, h);
    else return OUT_OF_RANGE;
}

// nearest neighbor interpolation
double DepthMap::nearest(const Vector2d pt, const int h) const
{
    int xd = xConv(pt[0]);
    int yd = yConv(pt[1]);
    if (isValid(xd, yd, h)) return at(xd, yd, h);
    else return OUT_OF_RANGE;
}


// nearest neighbor interpolation
double DepthMap::nearestSigma(const int u, const int v, const int h) const
{
    int xd = xConv(u);
    int yd = yConv(v);
    if (isValid(xd, yd, h)) return sigma(xd, yd, h);
    else return OUT_OF_RANGE;
}

// nearest neighbor interpolation
double DepthMap::nearestSigma(const Vector2d pt, const int h) const
{
    int xd = xConv(pt[0]);
    int yd = yConv(pt[1]);
    if (isValid(xd, yd, h)) return sigma(xd, yd, h);
    else return OUT_OF_RANGE;
}

// nearest neighbor interpolation
double DepthMap::nearestCost(const int u, const int v, const int h) const
{
    int xd = xConv(u);
    int yd = yConv(v);
    if (isValid(xd, yd, h)) return cost(xd, yd, h);
    else return OUT_OF_RANGE;
}

// nearest neighbor interpolation
double DepthMap::nearestCost(const Vector2d pt, const int h) const
{
    int xd = xConv(pt[0]);
    int yd = yConv(pt[1]);
    if (isValid(xd, yd, h)) return cost(xd, yd, h);
    else return OUT_OF_RANGE;
}


// to access the elements directly
double & DepthMap::at(const int x, const int y, const int h)
{
    return valVec[x + y*xMax + h*hStep];
}
const double & DepthMap::at(const int x, const int y, const int h) const
{
    return valVec[x + y*xMax + h*hStep];
}

// to access the elements directly
double & DepthMap::at(const int idx)
{
    return valVec[idx];
}
const double & DepthMap::at(const int idx) const
{
    return valVec[idx];
}

// to access the uncertainty directly
double & DepthMap::sigma(const int x, const int y, const int h)
{
    return sigmaVec[x + y*xMax + h*hStep];
}
const double & DepthMap::sigma(const int x, const int y, const int h) const
{
    return sigmaVec[x + y*xMax + h*hStep];
}

// to access the uncertainty directly
double & DepthMap::sigma(const int idx)
{
    return sigmaVec[idx];
}
const double & DepthMap::sigma(const int idx) const
{
    return sigmaVec[idx];
}

// to access the hypothesis cost directly
double & DepthMap::cost(const int x, const int y, const int h)
{
    return costVec[x + y*xMax + h*hStep];
}
const double & DepthMap::cost(const int x, const int y, const int h) const
{
    return costVec[x + y*xMax + h*hStep];
}

// to access the hypothesis cost directly
double & DepthMap::cost(const int idx)
{
    return costVec[idx];
}
const double & DepthMap::cost(const int idx) const
{
    return costVec[idx];
}

Vector2dVec DepthMap::getPointVec(const std::vector<int> & idxVec) const
{
    Vector2dVec result;
    result.reserve(idxVec.size());
    for (auto & idx : idxVec)
    {
        const int idxh = idx % hStep; //FIXME
        result.emplace_back(uConv(idxh % xMax), vConv(idxh / xMax));
    }
    return result;
}

Vector2dVec DepthMap::getPointVec() const
{
    Vector2dVec result;
    result.reserve(xMax * yMax);
    for (int y = 0; y < yMax; y++)
    {
        for (int x = 0; x < xMax; x++)
        {
            result.emplace_back(uConv(x), vConv(y)); 
        }
    }
    return result;
}

//TODO - Depecrated
void DepthMap::reconstructUncertainty(vector<int> & idxVec, 
            Vector3dVec & minDistVec, Vector3dVec & maxDistVec) const
{
    minDistVec.clear();
    maxDistVec.clear();
    idxVec.clear();
    vector<double> minVec;
    vector<double> maxVec;
    vector<int> idxBrutVec;
    for (int i = 0; i < valVec.size(); i++)
    {
        double d = valVec[i];
        if (d >= MIN_DEPTH)
        {
            double s = sigmaVec[i];
            // take d +- 2*sigma
            minVec.push_back(max(MIN_DEPTH, d - 2*s));
            maxVec.push_back(d + 2*s);
            idxBrutVec.push_back(i);
        }
    }
    
    Vector2dVec pointBrutVec = getPointVec(idxBrutVec);
    
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(pointBrutVec, reconstBrutVec, maskVec);
    
    for (int i = 0; i < reconstBrutVec.size(); i++)
    {
        if (maskVec[i])
        {
            Vector3d X = reconstBrutVec[i].normalized();
            minDistVec.push_back(X*minVec[i]);
            maxDistVec.push_back(X*maxVec[i]);
            idxVec.push_back(idxBrutVec[i]);        
        }
    }
}

//TODO - Depecrated
void DepthMap::reconstruct(vector<int> & idxVec, Vector3dVec & result) const
{
    result.clear();
    idxVec.clear();
    vector<double> depthVec;
    vector<int> idxBrutVec;
    for (int i = 0; i < valVec.size(); i++)
    {
        double d = valVec[i];
        if (d >= MIN_DEPTH)
        {
            depthVec.push_back(d);
            idxBrutVec.push_back(i);
        }
    }
    Vector2dVec pointBrutVec = getPointVec(idxBrutVec);
    
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(pointBrutVec, reconstBrutVec, maskVec);
    
    for (int i = 0; i < reconstBrutVec.size(); i++)
    {
        if (maskVec[i])
        {
            Vector3d X = reconstBrutVec[i].normalized();
            result.push_back(X*depthVec[i]);
            idxVec.push_back(idxBrutVec[i]);        
        }
    }
}

//TODO - Depecrated
void DepthMap::reconstruct(const Vector2dVec & queryPointVec,
        vector<int> & idxVec, Vector3dVec & result) const
{
    result.clear();
    idxVec.clear();
    Vector3dVec reconstBrutVec;
    vector<bool> maskVec;
    cameraPtr->reconstructPointCloud(queryPointVec, reconstBrutVec, maskVec);
    for (int i = 0; i < queryPointVec.size(); i++)
    {
        if (maskVec[i])
        {
            double d = nearest(queryPointVec[i]);
            if (d < MIN_DEPTH) continue;
            Vector3d X = reconstBrutVec[i].normalized();
            result.push_back(X*d);
            idxVec.push_back(i);        
        }
    }
}


vector<int> DepthMap::getIdxVec(const Vector2dVec & queryPointVec) const
{
    vector<int> outIdxVec;
    if(queryPointVec.size()!=0)
    {
        for(int i = 0; i < queryPointVec.size(); i++)
        {
            const int xd = xConv(queryPointVec[i][0]);
            const int yd = yConv(queryPointVec[i][1]);
            if (isValid(xd, yd)) outIdxVec.push_back( xd + yd*xMax );
            else outIdxVec.push_back(-1); // Idx of -1 to correspond to error
        }
    }
    else
    {
        for(int i=0; i < hStep; i++) // hStep === xMax * yMax
        {
            if (valVec[i] >= MIN_DEPTH) outIdxVec.push_back(i);
        }
    }
    return outIdxVec;
}


void DepthMap::reconstruct(MHPack & result, const uint32_t reconstFlags) const
{
    assert(not (reconstFlags & IMAGE_VALUES)); //FIXME imageValues are not implemented
    const int numHyps = (reconstFlags & ALL_HYPOTHESES) ? hMax : 1;

    // Convert query points to query indices
    vector<int> queryIdxVec;
    Vector2dVec queryPointVec;
    if (reconstFlags & QUERY_INDICES)
    {
        swap(queryIdxVec, result.idxVec);
        queryPointVec = getPointVec(queryIdxVec);
    }
    else if (reconstFlags & QUERY_POINTS)
    {
        queryIdxVec = getIdxVec(result.imagePointVec);
        swap(queryPointVec, result.imagePointVec);
    }
    else
    {
        queryIdxVec = getIdxVec();
        queryPointVec = getPointVec(queryIdxVec);
    }

    result.idxVec.clear();
    result.idxMapVec.clear();
    result.hypIdxVec.clear();
    result.costVec.clear();
    result.sigmaVec.clear();
    result.cloud.clear();
    result.valVec.clear();
//result.datatype = (minmax_flag) ? MHPack::MINMAX_DISTANCE_VEC_WITH_SIGN : MHPack::RECONSTRUCTION_WITH_SIGMA;

    vector<double> depthVec;
    for (int i = 0; i < queryIdxVec.size(); i++)
    {
        if (queryIdxVec[i] < 0 or queryIdxVec[i] > hStep) continue;

        for (int h = 0; h < numHyps; h++)
        {
            const int queryIdx = queryIdxVec[i];
            double depth = valVec[queryIdx + h*hStep];
            double sigma = sigmaVec[queryIdx + h*hStep]; //TODO discard points with sigma > sigmaMax
            if (depth < MIN_DEPTH)
            {
                if (reconstFlags & DEFAULT_VALUES and h == 0)
                {
                    depth = DEFAULT_DEPTH;
                    sigma = DEFAULT_SIGMA_DEPTH;
                }
                else continue;
            }
            
            if (reconstFlags & MINMAX) 
            {
                depthVec.push_back(max(depth - 3*sigma, MIN_DEPTH));
                depthVec.push_back(depth + 3*sigma);
            }
            else
            {
                depthVec.push_back(depth);
            }
            
            if (reconstFlags & SIGMA_VALUE) result.sigmaVec.push_back(sigma);
            result.idxVec.push_back(queryIdx);
            result.hypIdxVec.push_back(h);
            result.imagePointVec.push_back(queryPointVec[i]);
            if (reconstFlags & INDEX_MAPPING) result.idxMapVec.push_back(i);
        }
    }

    vector<bool> maskVec;
    Vector3dVec cloud;
    cameraPtr->reconstructPointCloud( result.imagePointVec, cloud, maskVec );
    
    result.cloud.reserve((reconstFlags & MINMAX) ? 2 * cloud.size() : cloud.size());
    for (int i = 0; i < cloud.size(); i++)
    {
        if (maskVec[i])
        {
            if (reconstFlags & MINMAX)
            {
                result.cloud.push_back(cloud[i].normalized() * depthVec[2*i]);
                result.cloud.push_back(cloud[i].normalized() * depthVec[2*i + 1]);
            }
            else
            {
                result.cloud.push_back(cloud[i].normalized() * depthVec[i]);
            }
        }
        else
        {
            //FIXME make all data in MHPack valid
            //that is discard all non-reconstructed points and other data associated to them
            if (reconstFlags & MINMAX)
            {
                result.cloud.emplace_back(0, 0, 0);
                result.cloud.emplace_back(0, 0, 0);
            }
            else
            {
                result.cloud.emplace_back(0, 0, 0);
            }
        }
    }
}


void DepthMap::project(const Vector3dVec & pointVec, Vector2dVec & result) const
{
    cameraPtr->projectPointCloud(pointVec, result);
}

void DepthMap::toMat(Mat32f & out) const
{
    out.create(yMax, xMax);
    copy(valVec.begin(), valVec.begin() + hStep, (float*)out.data);
}

void DepthMap::toInverseMat(Mat32f & out, const int layer) const
{
    out.create(yMax, xMax);
    const double* pInData = &valVec[0 + layer*hStep];
    float* pOutData = (float*)out.data;
    for (int i = 0; i < hStep; ++i)
    {
        if ( *pInData < 1e-3 ) *pOutData = 0;
        else *pOutData = 1 / (*pInData);
        ++pOutData;
        ++pInData;
    }
}

void DepthMap::sigmaToMat(Mat32f & out) const
{
    out.create(yMax, xMax);
    copy(sigmaVec.begin(), sigmaVec.begin() + hStep, (float*)out.data);
}

//TODO do not reconstruct all the points but a selected subset
// to avoid reconstruction of points with bad disparity
//TODO - Add support to insert multiple hypotheses into output depthmap
void DepthMap::wrapDepth(const DepthMap& dMap1, const DepthMap& dMap2,
        const Transformation<double> T12, DepthMap& output) const
{
	//Step 1 : Get point-cloud of first camera in first frame
	// vector<int> idx0Vec;
	// Vector3dVec cloud11;
    MHPack cloud11MH;
	// dMap1.reconstruct(idx0Vec, cloud11);
    dMap1.reconstruct(cloud11MH);

	//Step 2 : Transform above into second frame
	Vector3dVec cloud12;
	T12.inverseTransform(cloud11MH.cloud, cloud12);

	//Step 3 : Reproject points into second camera
	Vector2dVec point12Vec;
	dMap2.project(cloud12, point12Vec);

	//Step 4 : For reprojected points, reconstruct point-cloud of second camera in second frame
	// Vector3dVec cloud22;
	// vector<int> idx1Vec;
    MHPack cloud22MH;
	// dMap2.reconstruct(point12Vec, idx1Vec, cloud22);
	// Do not reconstruct all hypotheses, we need a 1-1 map between input and output
    dMap2.reconstruct(cloud22MH, QUERY_POINTS); 

	//Step 5 : Transform above into first frame
	Vector3dVec cloud21;
	T12.transform(cloud22MH.cloud, cloud21);

	//Step 6 : Project above points along corresponding depth vectors
    output = dMap1;
    output.setTo(0, 1);
	// for(int i = 0; i < idx1Vec.size(); ++i)
    for (int i = 0; i < cloud22MH.idxVec.size(); i++)
	{
	    // int idx1 = idx1Vec[i];
	    // int idx0 = idx0Vec[idx1];
        const int idx1 = cloud22MH.idxVec[i];
        const int idx0 = cloud11MH.idxVec[idx1];
		const Vector3d & X2 = cloud21[i];
		// const Vector3d & X1 = cloud11[idx1];
        const Vector3d & X1 = cloud11MH.cloud[idx0];
		//Calculate dot-product to get the distance as the projection along the line
		output.at(idx0) =  X2.dot(X1.normalized());
		// output.sigma(idx0) = dMap2.nearestSigma(point12Vec[idx1]);
        output.sigma(idx0) = cloud22MH.valVec[idx1]; // As we used RECONSTRUCTION_WITH_SIGMA
	}
}

//TODO - Add support to insert multiple hypotheses into output depthmap
DepthMap DepthMap::wrapDepth(const Transformation<double> T12) const
{
    DepthMap dMap2(cameraPtr, *this);

    //Step 1 : Get point-cloud of current frame
    MHPack cloud11MH;
    reconstruct(cloud11MH, 0);

    //Step 2 : Transform cloud into keyframe
    Vector3dVec cloud12;
    T12.inverseTransform(cloud11MH.cloud, cloud12);

    //Step 3 : Project point cloud on keyframe
    Vector2dVec point12Vec;
    cameraPtr->projectPointCloud(cloud12, point12Vec);

    //Step 4 : Fill in data for depthmap
    vector<int> idx12Vec = getIdxVec(point12Vec);
    for (int i = 0; i < idx12Vec.size(); i++)
    {
        const int idx2 = idx12Vec[i];
        const int idx1 = cloud11MH.idxVec[i];
        if(idx2 != -1)
        {
            dMap2.at(idx2) = cloud12[i].norm();
            dMap2.sigma(idx2) = sigma(idx1);
            dMap2.cost(idx2) = cost(idx1);
        }
    }

    return dMap2;
}

DepthMap DepthMap::generatePlane(const ICamera * camera, const ScaleParameters & params, 
        Transformation<double> TcameraPlane, const Vector3dVec & polygonVec)
{
    DepthMap depth(camera, params);
    Vector3d t = TcameraPlane.trans();
    Vector3d z = TcameraPlane.rotMat().col(2);
    Vector3dVec polygonCamVec;
    TcameraPlane.transform(polygonVec, polygonCamVec);
    for (int v = 0; v < params.yMax; v++)
    {
        for (int u = 0; u < params.xMax; u++)
        {
            depth.at(u, v) = OUT_OF_RANGE;
            depth.sigma(u, v) = OUT_OF_RANGE;
            Vector3d vec; // the direction vector
            if (not camera->reconstructPoint(Vector2d(params.uConv(u), params.vConv(v)), vec)) continue;
            double zvec = z.dot(vec);
            if (zvec < 1e-3) 
            {
                continue;
            }
            bool inside = true;
            for (int i = 0; i < polygonCamVec.size(); i++)
            {
                int j = (i + 1) % polygonCamVec.size();
                Vector3d normal = polygonCamVec[i].cross(polygonCamVec[j]);
                if (vec.dot(normal) < 0)
                {
                    inside = false;
                    break;
                }
            }
            if (not inside) continue;
            double tz = t.dot(z);
            double alpha = tz / zvec;
            vec *= alpha;
            depth.at(u, v) = vec.norm();
            depth.sigma(u, v) = 1;
        }
    }
    return depth;
}

void DepthMap::pixelMedianFilter(const int x, const int y, const int h, DepthMap & dst)
{
    vector<double> depths;
    vector<double> sigmas;
    for (int y1 = max(0, y-1); y1 < min(yMax, y+1); ++y1)
    {
        for (int x1 = max(0, x-1); x1 < min(xMax, x+1); ++x1)
        {
            for (int h1 = 0; h1 < hMax; ++h1)
            {
                const double d = at(x1, y1, h1);
                if (d > MIN_DEPTH)
                {
                    depths.push_back(d);
                    sigmas.push_back(sigma(x1, y1, h1));
                }
            }
        }
    }
    if (depths.size() < 1) return; // Do nothing
    // Sort vectors in order of increasing depth
    for (int i = 0; i < depths.size()-1; ++i)
    {
        for (int j = i+1; j < depths.size(); ++j)
        {
            if(depths[i] > depths[j])
            {
                const double temp_d = depths[i];
                depths[i] = depths[j];
                depths[j] = temp_d;
                const double temp_s = sigmas[i];
                sigmas[i] = sigmas[j];
                sigmas[j] = temp_d;
            }
        }
    }
    //Choose median element and assign it to point
    const int median = depths.size() / 2;
    dst.at(x, y, h) = depths[median];
    dst.sigma(x, y, h) = depths[median];
}

void DepthMap::pixelAverageFilter(const Vector3iVec & matches, DepthMap & dst)
{
    if (matches.size() == 0)
    {
        std::cerr << "ERROR: Attempted to use pixelAverageFilter with zero matches." << std::endl;
        return; // return to parent function, don't crash program
    }
    double sum_d = 0, sum_s = 0;
    for (int i = 0; i < matches.size(); ++i)
    {
        const Vector3i & point = matches[i];
        sum_d += at(point[0], point[1], point[2]);
        sum_s += sigma(point[0], point[1], point[2]);
    }
    const Vector3i & base = matches[0];
    dst.at(base[0], base[1], base[2]) = sum_d / matches.size();
    dst.sigma(base[0], base[1], base[2]) = sum_s / matches.size();
}

//works only for the first hypothesis so far
//TODO implement for multi-hyp case
void DepthMap::filterNoise()
{
    const int minMatches = 2;

    // Create copy of current depthmap, to avoid flow
    DepthMap myCopy = *this;

    const array<int, 8> dxArr =  {-1,  0,  1,  1,  1,  0, -1, -1};
    const array<int, 8> dyArr  = { 1,  1,  1,  0, -1, -1, -1,  0};
    
    const int CENTRAL_WEIGHT = 5;
    // Three loops to loop through every hypothesis
    for (int y = 1; y < yMax - 1; ++y)
    {
        for (int x = 1; x < xMax - 1; ++x)
        {
            double depthVal = myCopy.at(x, y);
            if (depthVal == OUT_OF_RANGE) continue;
            double sigmaVal = myCopy.sigma(x, y);
            int countFilled = 0, countMatches = 0;
            double acc = depthVal * CENTRAL_WEIGHT;
            for (int i = 0; i < 8; i++)
            {
                double neighDepthVal = myCopy.at(x + dxArr[i], y + dyArr[i]);
//                double neighSigmaVal = myCopy.sigma(x + dxArr[i], y + dyArr[i]);
                if (neighDepthVal == OUT_OF_RANGE) continue;
                countFilled++;
                if (abs(depthVal - neighDepthVal) > 2*sigmaVal) continue;
                countMatches++;
                acc += neighDepthVal;
            }
            if (countMatches < minMatches and countMatches < countFilled or countFilled < 2)
            {
                at(x, y) = OUT_OF_RANGE;
                sigma(x, y) = OUT_OF_RANGE;
            }
            else
            {
                at(x, y) = acc / (countMatches + CENTRAL_WEIGHT);
            }
        }
    }

}

void DepthMap::merge(const DepthMap & depth2)
{
    assert((ScaleParameters)(*this) == (ScaleParameters)depth2);

    // Filter merge all new hypotheses
    for (int y = 0; y < yMax; y++)
    {
        for (int x = 0; x < xMax; x++)
        {
            for (int h2 = 0; h2 < depth2.hMax; h2++)
            {
                const double d2 = depth2.at(x, y, h2);
                if (d2 < MIN_DEPTH) continue;
                const double sigma2 = depth2.sigma(x, y, h2);
                if (sigma2 > 3) continue;
                filterPushHypothesis(x, y, d2, sigma2);
            }
        }
    }

    // Remove bad hypotheses
    costRejection();
}

void DepthMap::regularize()
{
    const int minMatches = 2;
    // Three loops to loop through every hypothesis
    for (int h = 0; h < hMax - 1; ++h)
    {
        for (int y = 0; y < yMax; ++y)
        {
            for (int x = 0; x < xMax; ++x)
            {
                if (at(x, y, h) < MIN_DEPTH)
                {
                    int h2 = h + 1;
                    while (h2 < hMax and at(x, y, h2) < MIN_DEPTH) h2++;
                    if (h2 == hMax) continue;
                    at(x, y, h) = at(x, y, h2);
                    sigma(x, y, h) = sigma(x, y, h2);
                    cost(x, y, h) = cost(x, y, h2);
                    at(x, y, h2) = OUT_OF_RANGE;
                }             
            }
        }
    }
}
