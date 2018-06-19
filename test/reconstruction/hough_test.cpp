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

#include "io.h"
#include "ocv.h"
#include "eigen.h"
#include "json.h"

#include "geometry/geometry.h"
#include "utils/curve_rasterizer.h"
#include "projection/eucm.h"





Polynomial2 computePolynomial(Vector3d plane, const double * const params)
{
    Polynomial2 surf;
    
    const double alpha = params[0];
    const double beta = params[1];
    const double fu = params[2];
    const double fv = params[3];
    const double u0 = params[4];
    const double v0 = params[5];
    
    // intermediate variables
    const double gamma = 1 - alpha;
    const double ag = (alpha - gamma);
    const double a2b = alpha*alpha*beta;
    const double fufv = fu * fv;
    const double fufu = fu * fu;
    const double fvfv = fv * fv;
    const double & A = plane[0];
    const double & B = plane[1];
    const double & C = plane[2];
    const double AA = A * A;
    const double BB = B * B;
    const double CC = C * C;
    const double CCfufv = CC * fufv;
    const double dd = CCfufv / (AA + BB); //squared distance to the projection center in pixels
    if ((AA + BB) > 0 and dd < 1.) // the curve passes through the projection center
    {
        
        surf.kuu = surf.kuv = surf.kvv = 0;
        surf.ku = A / fu;
        surf.kv = B / fv;
        
        const double normABinv = 1./sqrt(AA + BB);
        const double Cnorm = C / sqrt(AA + BB + CC);
        const double du = -A * Cnorm * normABinv * fu;
        const double dv = -B * Cnorm * normABinv * fv;
        surf.k1 = -(u0  + du)* A/ fu - (v0 + dv) * B / fv;
        
//        surf.k1 = -u0 * A/ fu - v0 * B / fv;
    }
    else
    {
        // compute first 4 coefficients directly
        surf.kuu = (AA*ag + CC*a2b)/(CC*fufu);  // kuu
        surf.kuv = 2*A*B*ag/(CCfufv);  // kuv
        surf.kvv = (BB*ag + CC*a2b)/(CC*fvfv);  // kvv
        surf.ku = 2*(-(AA*fv*u0 + A*B*fu*v0)*ag - 
                        A*C*fufv*gamma - CC*a2b*fv*u0)/(CCfufv*fu);  // kv
        surf.kv = 2*(-(BB*fu*v0 + A*B*fv*u0)*ag - 
                        B*C*fufv*gamma - CC*a2b*fu*v0)/(CCfufv*fv);  // kv
                        
        // the last one is computed using the fact that
        // the epipolar curves pass through the epipole
        surf.k1 = (A*A*fv*fv*u0*u0*(alpha - gamma) + 2*A*B*fu*fv*u0*v0*(alpha - gamma) 
                    + 2*A*C*fu*fv*fv*gamma*u0 + B*B*fu*fu*v0*v0*(alpha - gamma) + 2*B*C*fu*fu*fv*gamma*v0 
                    + C*C*alpha*alpha*beta*fu*fu*v0*v0 + C*C*alpha*alpha*beta*fv*fv*u0*u0 - C*C*fu*fu*fv*fv)/(C*C*fu*fu*fv*fv);

    }
    return surf;
}

const int MARGIN = 5;
const double GRAD_THRESH = 0.001;
const int SEARCH_SIZE = 2;
const double ACC_THRESH = 7;
const double ACC_DELTA_THRESH = 0.05;

// solves equation z = 0
void findCircleIntersections(const EnhancedCamera & camera, const Vector3d & n, Vector2d & p1, Vector2d & p2)
{
    double r = 1. / ( camera.getAlpha() * sqrt(camera.getBeta()) );
    Vector2d pn(n[0], n[1]);
    if (abs(n[0]) + abs(n[1]) < 1e-4) 
    {
        p1[0] = r;
        p1[1] = 0;
        
        p2[0] = 0;
        p2[1] = r;
        return;
    }
    
    pn.normalize();
    pn *= r;
    
    p1[0] = pn[1] * camera.getFocalU() + camera.getCenterU();
    p1[1] = -pn[0] * camera.getFocalV() + camera.getCenterV();
    
    p2[0] = -pn[1] * camera.getFocalU() + camera.getCenterU();
    p2[1] = pn[0] * camera.getFocalV() + camera.getCenterV();
}

//TODO write a class




bool solveQuadratic(double A, double B, double C, double & x1, double & x2)
{
    if (abs(A) < 1e-8)
    {
        x1 = -C / B;
        x2 = -DOUBLE_MAX;
        return true;
    }
    double Delta = B * B - 4 * A * C;
    cout << "DELTA : "<< Delta << endl;
    if (Delta < 0) return false;
    double Ainv = 0.5 / A;
    double sqrtDelta = sqrt(Delta);
    x1 = (-B + sqrtDelta) * Ainv;
    x2 = (-B - sqrtDelta) * Ainv;
    return true;
}

//fixed U coordinate
bool solvePolyU(const Polynomial2 & poly, const double u, double & v1, double & v2)
{
    double B = poly.kuv * u + poly.kv;
    double C = poly.kuu * u * u + poly.ku * u + poly.k1;
    cout << poly.kvv << "   " << B << "   " << C << endl;
    bool res = solveQuadratic(poly.kvv, B, C, v1, v2);
    return res;
}

//fixed V coordinate
bool solvePolyV(const Polynomial2 & poly, const double v, double & u1, double & u2)
{
    double B = poly.kuv * v + poly.ku;
    double C = poly.kvv * v * v + poly.kv * v + poly.k1;
    cout << poly.kuu << "   " << B << "   " << C << endl;
    bool res = solveQuadratic(poly.kuu, B, C, u1, u2);
    return res;
}

const double NOT_VALID = -1;
const double BOTH_VALID = -2;
//among the two given values keeps the one which 
double chooseValidSolution(const double & valMax, double x1, double x2)
{
    x1 = round(x1);
    x2 = round(x2);
    bool valid1 = (x1 >= 0 and x1 < valMax);
    bool valid2 = (x2 >= 0 and x2 < valMax); 
    if (valid1 and valid2)
    {
        return BOTH_VALID;
    }
    else if (valid1) return x1;
    else if (valid2) return x2;
    else return NOT_VALID;
}

Vector2d checkPolySolution(const Polynomial2 & poly, const Vector2d & pt, const Vector2d & pt1, const Vector2d & pt2)
{
    double val1 = abs( poly(0.5*(pt[0] + pt1[0]), 0.5*(pt[1] + pt1[1])) );
    double val2 = abs( poly(0.5*(pt[0] + pt2[0]), 0.5*(pt[1] + pt2[1])) );
    cout << "VALS1 : " << val1 << "  VAL2 : " << val2 << endl;
    cout << poly(pt[0], pt[1]) << endl;
    cout << poly(pt1[0], pt1[1]) << endl;
    cout << poly(pt2[0], pt2[1]) << endl;
    if (val1 < val2) return pt1;
    else return pt2;
}

double midpointCheck(const double & x0, const double & x1, const double & x2)
{
    if (abs(x1 - x0) < abs(x2 - x0)) return x1;
    else return x2;
}

bool distanceCheck(const EnhancedCamera & camera, Vector2d pt)
{
    double r2 = camera.getFocalU() * camera.getFocalV() / 
                ( camera.getAlpha() * camera.getAlpha() * camera.getBeta() );
    pt[0] -= camera.getCenterU();
    pt[1] -= camera.getCenterV();
    return (pt.squaredNorm() < r2);
}

Vector2d findPoint(const EnhancedCamera & camera, const Polynomial2 & poly, const Vector2d & pt, const double u0, const double v0)
{
    double v1, v2;
    Vector2d candidate1;
    if (solvePolyU(poly, u0, v1, v2))
    {
        cout << "v1 : " << v1 << "   v2 : " << v2 << endl; 
        double res = chooseValidSolution(camera.height, v1, v2);
        if (res >= 0) candidate1 = Vector2d(u0, res);
        else if (res == BOTH_VALID)
        {
            candidate1 = Vector2d(u0, midpointCheck(camera.getCenterV(), v1, v2));
            
        }
        if (distanceCheck(camera, candidate1)) return candidate1;
    }
    double u1, u2;
    Vector2d candidate2;
    if (solvePolyV(poly, v0, u1, u2))
    {
        cout << "U1 : " << u1 << "   U2 : " << u2 << endl; 
        double res = chooseValidSolution(camera.width, u1, u2);
        if (res >= 0) candidate2 = Vector2d(res, v0);
        else if (res == BOTH_VALID)
        {
            candidate2 = Vector2d(midpointCheck(camera.getCenterU(), u1, u2), v0);
        }
        if (distanceCheck(camera, candidate2)) return candidate2;
    }
    throw;
}

bool bringToBorder(const EnhancedCamera & camera, const Polynomial2 & poly, Vector2d & pt)
{
    pt[0] = round(pt[0]);
    pt[1] = round(pt[1]);
    
    if (pt[0] < camera.width    and pt[0] >= 0 and
        pt[1] < camera.height   and pt[1] >= 0) return true;
        
    if (pt[0] < camera.width / 2 and pt[1] < camera.height / 2)
    {
        pt = findPoint(camera, poly, pt, 0, 0);
        return true;
    }
    else if (pt[0] < camera.width / 2)
    {
        pt = findPoint(camera, poly, pt, 0, camera.height - 1);
        return true;
    }
    else if (pt[1] < camera.height / 2)
    {
        pt = findPoint(camera, poly, pt, camera.width - 1, 0);
        return true;
    }
    else
    {
        pt = findPoint(camera, poly, pt, camera.width - 1, camera.height - 1);
        return true;
    }
}

//camera must contain Width/height information
bool findTerminalPoints(const EnhancedCamera & camera, const Vector3d & norm, const Polynomial2 & poly, Vector2d & pt1, Vector2d & pt2)
{
	//find circle intersections
    findCircleIntersections(camera, norm, pt1, pt2);
    cout << pt1.transpose() << "   " << pt2.transpose() << endl;
    bringToBorder(camera, poly, pt1);
    bringToBorder(camera, poly, pt2);
    cout << pt1.transpose() << "   " << pt2.transpose() << endl;
}


Vector2d subpixel(const Mat32f & acc, int u, int v)
{
//    Matrix<double, 8, 5> A;
//    Matrix<double, 8, 1> B;
    double  VAL_THRESH = acc(v, u);
    double uAcc = 0, vAcc = 0, numAcc = 0;
    const double AVG_SIZE = 1;
    for (int du = -AVG_SIZE ; du <= AVG_SIZE ; du++)
    {
        for (int dv = -AVG_SIZE ; dv <= AVG_SIZE ; dv++)
        {   
            double val = acc(v + dv, u + du);
            numAcc += val;
            uAcc += du * val;
            vAcc += dv * val;
//            if (du == 0 and dv == 0) continue;
//            A(i, 0) = du;
//            A(i, 1) = dv;
//            A(i, 2) = du * du;
//            A(i, 3) = du * dv;
//            A(i, 4) = dv * dv;
//            B(i) = acc(v + dv, u + du) - k1;
        }
    }
    return Vector2d(u + uAcc / numAcc, v + vAcc / numAcc);
}

//Vector2d subpixel(const Mat32f & acc, int u, int v)
//{
////    Matrix<double, 8, 5> A;
////    Matrix<double, 8, 1> B;
//    double  k1 = acc(v, u);
//    double uAcc = 0, vAcc = 0, numAcc = 0;
//    const double AVG_SIZE = 1;
//    for (int du = -AVG_SIZE ; du <= AVG_SIZE ; du++)
//    {
//        for (int dv = -AVG_SIZE ; dv <= AVG_SIZE ; dv++)
//        {   
//            double val = acc(v + dv, u + du);
//            numAcc += val;
//            uAcc += du * val;
//            vAcc += dv * val;
////            if (du == 0 and dv == 0) continue;
////            A(i, 0) = du;
////            A(i, 1) = dv;
////            A(i, 2) = du * du;
////            A(i, 3) = du * dv;
////            A(i, 4) = dv * dv;
////            B(i) = acc(v + dv, u + du) - k1;
//        }
//    }
//    return Vector2d(u + uAcc / numAcc, v + vAcc / numAcc);
//}

int main(int argc, char** argv) 
{

    ptree root;
    read_json(argv[1], root);
    
    EnhancedCamera camera(  root.get<int>("width"), 
                            root.get<int>("height"),
                            readVector<double>(root.get_child("camera_params")).data() );
                            
    EnhancedCamera accCam( readVector<double>(root.get_child("acc_params")).data() );
    
    vector<string> fnameVec;
    for (auto & x : root.get_child("names"))
    {
        fnameVec.emplace_back(x.second.get_value<string>());
    }
    int count = 0;
    for (auto & name : fnameVec)
    {
        count++;
        Mat8u img = imread(name, 0);
        Mat32f gx, gy, gnorm;
        Sobel(img, gx, CV_32F, 1, 0, 3, 1./2048);
        Sobel(img, gy, CV_32F, 0, 1, 3, 1./2048);
//        Sobel(img, gx, CV_32F, 1, 0, 1, 1./512);
//        Sobel(img, gy, CV_32F, 0, 1, 1, 1./512);
        gnorm = gx.mul(gx) + gy.mul(gy);
        
        Mat32f acc(2 * accCam.getCenterV(), 2 * accCam.getCenterU());
        acc.setTo(0);
        
        
        const double alpha = camera.getParams()[0];
        const double beta = camera.getParams()[1];
        const double gamma = 1. - alpha;
        
        imshow("grad", gnorm*5);
//        waitKey();
        
        for (int v = MARGIN; v < gnorm.rows - MARGIN; v++)
        {
            for (int u = MARGIN; u < gnorm.cols - MARGIN; u++)
            {
                if (gnorm(v, u) < GRAD_THRESH) continue;
                
                Vector3d X;
                camera.reconstructPoint(Vector2d(u, v), X);
                const double & x = X[0];
                const double & y = X[1];
                const double & z = X[2];
                
                double Kx = 2 * alpha * alpha * beta * x;
                double Ky = 2 * alpha * alpha * beta * y;
                
                double P = 2 * (alpha - gamma);
                double Q = 2 * gamma;
                
                double fx = gx(v, u);
                double fy = gy(v, u);
                
                Vector3d Y( -fy * (P * z + Q),
                            fx * (P * z + Q),
                            fy * Kx - fx * Ky);
                
                Vector3d ABC = X.cross(Y);
                if (ABC[2] < 0) ABC *= -1;
                Vector2d abc;
                accCam.projectPoint(ABC, abc);
                
    //            cout    << v << "   " << u << "  #   " << X.transpose() << "   #    "  
    //                    << ABC.transpose() << "   #   " << abc.transpose() << endl;
                
                acc(round(abc[1]), round(abc[0]))++;
                
                
                if ( ABC[2] * ABC[2] / (ABC[0] * ABC[0] + ABC[1] * ABC[1]) < 3e-3)
                {
                    accCam.projectPoint(-ABC, abc);
                    acc(round(abc[1]), round(abc[0]))++;
                }
                
    //            CurveRasterizer<int, Polynomial2> raster1(Vector2i(u, v), Vector2i(u + 100*fy, v - 100*fx),
    //                             computePolynomial(ABC, camera.getParams()));
    //            
    //            CurveRasterizer<int, Polynomial2> raster2 = raster1;                
    //            Mat8uc3 imgTrace;
    //            cvtColor(img, imgTrace, CV_GRAY2BGR);
    //            cv::circle(imgTrace, Point(u, v), 3, Scalar(255, 100, 100), 2); 
    //            const int LENGTH = 150;
    //            for (int i = 0; i < LENGTH; i++)
    //            {
    //                cv::circle(imgTrace, Point(raster1.u, raster1.v), 1, Scalar(0, 255, 0), 2);
    //                raster1.step();
    //                cv::circle(imgTrace, Point(raster2.u, raster2.v), 1, Scalar(0, 255, 0), 2);
    //                raster2.unstep();
    //            }
    //    
    //            imshow("accumulator", imgTrace);
    //            waitKey();
            }
        }
        
        
        //Find local maxima
        GaussianBlur(acc, acc, Size(5, 5), 0.9, 0.9);
        GaussianBlur(gnorm, gnorm, Size(5, 5), 0.9, 0.9);
        imshow("accumulator", acc / 100);
//        waitKey();    
        
        Mat8uc3 imgTrace;
        cvtColor(img, imgTrace, CV_GRAY2BGR);
        double radius = camera.getFocalU() / ( camera.getAlpha() * sqrt(camera.getBeta()) );
                    cv::circle(imgTrace, Point(camera.getCenterU(), camera.getCenterV()), radius, Scalar(100, 100, 255), 2); 
                    
        Mat8uc3 accTrace(2 * accCam.getCenterV(), 2 * accCam.getCenterU());
        accTrace.setTo(0);
        
        for (int v = SEARCH_SIZE; v < acc.rows - SEARCH_SIZE; v++)
        {
            for (int u = SEARCH_SIZE; u < acc.cols - SEARCH_SIZE; u++)
            {
                double val = acc(v, u);
                if (val < ACC_THRESH) continue;
                
                bool isMax = true;
                for (int du = -SEARCH_SIZE; du <= SEARCH_SIZE and isMax; du++)
                {
                    for (int dv = -SEARCH_SIZE; dv <= SEARCH_SIZE and isMax; dv++)
                    {
                        if (du == 0 and dv == 0) continue;
                        
                        if (val <= acc(v + dv, u + du) + ACC_DELTA_THRESH) 
                        {
                            isMax = false;
                        }
                    }
                }
                if (isMax) //trace
                {
                    Vector3d ABC;
                    accCam.reconstructPoint(subpixel(acc, u, v), ABC);
                    if (ABC[2] < 0) continue;
                    cout << "DETECTED " << name << endl;
                    
//                    //Perpendicual Directions
//                    Vector3d d1(ABC[1] + ABC[2], -ABC[0], -ABC[0]);
//                    if (d1.squaredNorm() < 1e-3) d1 = Vector3d(ABC[1] - ABC[2], -ABC[0], ABC[0]);
//                    
//                    Vector3d d2(-ABC[1] , ABC[0] + ABC[2], -ABC[1]);
//                    if (d2.squaredNorm() < 1e-3) d2 = Vector3d(-ABC[1], ABC[0] - ABC[2], ABC[1]);
//                    
                    Vector2d pt1, pt2;
//                    
//                    camera.projectPoint(d1, p1);
//                    camera.projectPoint(d2, p2);
                    Polynomial2 poly = computePolynomial(ABC, camera.getParams());
                    
                    findTerminalPoints(camera, ABC, poly, pt1, pt2);
                    CurveRasterizer<int, Polynomial2> raster1(round(pt1), round(pt2), poly);
                    
//                    CurveRasterizer<int, Polynomial2> raster2 = raster1;                
                    
                    cv::circle(accTrace, Point(u, v), 3, Scalar(255, 100, 100), 2); 
                    
                    cv::circle(imgTrace, Point(pt1[0], pt1[1]), 4, Scalar(255, 100, 100), 2); 
                    cv::circle(imgTrace, Point(pt2[0], pt2[1]), 4, Scalar(255, 100, 100), 2); 
                    
                    
                    
                    const int LENGTH = 1500;
                    for (int i = 0; i < LENGTH and abs(raster1.u - pt2[0]) + abs(raster1.v - pt2[1]) > 2; i++)
                    {
//                        Vector3d Xrec;
//                        camera.reconstructPoint(Vector2d(raster1.u, raster1.v), Xrec);
//    //                    if ()
//                        if (Xrec.dot(ABC) < 1e-2
//                            and raster1.u >= 0 and raster1.u < gnorm.cols 
//                            and raster1.v >= 0 and raster1.v < gnorm.rows
//                            and gnorm(raster1.v, raster1.u) > 0.5*GRAD_THRESH)
//                        {
//                            
//                            cv::circle(imgTrace, Point(raster1.u, raster1.v), 0, Scalar(0, 255, 0), 1);
                            imgTrace(raster1.v, raster1.u) = cv::Vec3b(0, 255, 0);
//                        }
                        raster1.step();
//                        
//                        camera.reconstructPoint(Vector2d(raster2.u, raster2.v), Xrec);
//                        if (Xrec.dot(ABC) < 1e-2
//                            and raster2.u >= 0 and raster2.u < gnorm.cols 
//                            and raster2.v >= 0 and raster2.v < gnorm.rows
//                            and gnorm(raster2.v, raster2.u) > 0.5*GRAD_THRESH)
//                        {
//                            cv::circle(imgTrace, Point(raster2.u, raster2.v), 1, Scalar(0, 255, 0), 2);
//                        }
//                        raster2.unstep();
                    }
                    imshow("accTrace", accTrace);
                    imshow("res", imgTrace);
                    imwrite("res_" + to_string(count) + ".png", imgTrace);
                    waitKey();
                    
                }
            }
            
        }
        imshow("accTrace", accTrace);
        imshow("res", imgTrace);
        imwrite("res_" + to_string(count) + ".png", imgTrace);
        waitKey();
    }
}



