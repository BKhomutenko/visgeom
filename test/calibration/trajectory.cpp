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

#include "geometry/geometry.h"
#include "calibration/trajectory_generation.h"
#include "projection/eucm.h"

const int PARAM_NUM = 3;
class CircularTrajectory : public ITrajectory
{
public:
    CircularTrajectory(int numberSteps) : 
        _numberSteps(numberSteps) {}
        
    virtual void compute(const double * params, vector<Transf> & trajVec, 
            vector<Matrix6d> & covVec) const
    {
        trajVec.clear();
        
        //control covariance
        Matrix2d covVW;
        //TODO make parameters
        
        Matrix6d covAbs = Matrix6d::Identity()*1e-2;
        double alpha = params[0]; //turn angle
        double dist = 0.05;//params[PARAM_NUM*circle + 1];
        double ca = cos(alpha*0.5);
        double sa = sin(alpha*0.5);
        Transf xi0(params[2], 0, 0,
                     0, 0, params[1]);
        //TODO modelizer proprement l'odometrie            
        covVW <<    1e-4,        0, 
                    0,       1e-4;
                
        //circular motion model
        Transf dxi(dist * ca, dist * sa, 0, 0, 0, alpha);
        //motion jacobian
        Matrixd<6, 2> dxidu;
        dxidu <<    ca,     -dist/2 * sa,
                    sa,     dist/2 * ca,
                    0,          0,
                    0,          0,
                    0,          0,
                    0,          1;
        Matrix6d covIncr = dxidu * covVW * dxidu.transpose();
        trajVec.push_back(xi0);
        covVec.push_back(covAbs);
        
        //Screw transformation matrix
        Matrix6d L = dxi.screwTransfInv();
        Matrix6d covOdom = covIncr;
        for (int i = 1; i < _numberSteps; i++)
        {
            trajVec.push_back(trajVec.back().compose(dxi));
            covOdom = L * covOdom * L.transpose() + covIncr;
            covVec.push_back(covOdom + covAbs);
        }
    }
    
    virtual int paramSize() const { return  PARAM_NUM; }
    
    int _numberSteps;
};


int main(int argc, char** argv) 
{
    
    int circleCount = 2;
    vector<ITrajectory*> trajVec;
    double numberSteps = 30;
    vector<double> paramVec;
    for (int i = 0; i < circleCount; i++)
    {
        paramVec.push_back(i * 0.001 + 0.01);
//        paramVec.push_back(0.01);
//        paramVec.push_back(0);
//        paramVec.push_back(0);
        paramVec.push_back(-0.1);
        paramVec.push_back(1);
        trajVec.push_back(new CircularTrajectory(numberSteps));
    }
    
    vector<Transf> xiOdomVec;
    vector<Matrix6d> covOdomVec;
//    traj->compute(paramVec.data(), xiOdomVec, covOdomVec);
//    
//    for (int i = 0; i < xiOdomVec.size(); i++)
//    {
//        cout << xiOdomVec[i] << endl;
//        cout << covOdomVec[i] << endl << endl;
//    }
    
    
    
    //generate the board
    Vector3dVec board;
    const int N = 5, M = 8;
    double step = 0.1;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            board.emplace_back(0, step * j, step * i);
        }
    }
    
    //ZOE
    /*
    Transf xiBoard(7, -0.4, 0.5, 0, 0, 0);
    Matrix3d R;
    R <<    0,      0,      1, 
            -1,     0,      0,
            0,      -1,     0;
    Transf xiCam(Vector3d(3, 0.35, 0.3), R);
    cout << xiCam.rot() << endl;
    */
    
    //PIONEER
    Transf xiBoard(3, -0.4, 0.5, 0, 0, 0);
    Matrix3d R;
    R <<    0,      0,      1, 
            -1,     0,      0,
            0,      -1,     0;
    Transf xiCam(Vector3d(0.3, 0, 0.3), R);
    cout << xiCam.rot() << endl;
    
    
//    Transf xiCam2(0, 0, 0, 0, 0, 0);
    
    vector<double> camParams{0.719981,  1.03894,  381.974,  382.378,  523.6,  366.235};
    EnhancedCamera cam(1024, 768, camParams.data());
    TrajectoryVisualQuality * quality = new TrajectoryVisualQuality(
                                                trajVec, &cam, xiCam, xiBoard, board,
                                                Matrix6d::Identity(), Matrix2d::Identity());
    
    //////////////////////////////////
    //Test the visual covariance
    //////////////////////////////////
    Matrix6d C = quality->visualCov(Transf(0.1, 0.1, 0, 0, 0, .1).compose(xiCam));
    cout << C << endl << endl;
    cout << C.inverse() << endl;
    
    cout << quality->imageLimitsCost(Transf(0.1, 0.1, 0, 0, 0, 1.2).compose(xiCam)) << endl;
    cout << quality->imageLimitsCost(Transf(0, 0, 0, 0, 0, 0).compose(xiCam)) << endl;
    
    ceres::GradientProblem problem(quality);

    ceres::GradientProblemSolver::Options options;
    options.max_num_iterations = 1000;
    options.minimizer_progress_to_stdout = true;
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, paramVec.data(), &summary);
    
    std::cout << summary.FullReport() << endl;
    
    for (int i = 0; i < circleCount; i++)
    {
        for (int j = 0; j < PARAM_NUM; j++)
        {
             cout << paramVec[i * PARAM_NUM + j] << "   ";
        }
        
        cout <<  endl;
    }
    ofstream myfile;
    myfile.open ("/home/bogdan/projects/python/trajectory/traj3");
    
    
    
    for (int trajIdx = 0; trajIdx < trajVec.size(); trajIdx++)
    {
        auto traj = trajVec[trajIdx];
        traj->compute(paramVec.data() + trajIdx * traj->paramSize(), xiOdomVec, covOdomVec);
        
        for (int i = 0; i < xiOdomVec.size(); i++)
        {
            myfile << xiOdomVec[i].trans()[0] << "   " << xiOdomVec[i].trans()[1]  << endl;
            cout << xiOdomVec[i] << endl;
//            cout << quality->imageLimitsCost(xiOdomVec[i].compose(xiCam)) << endl;
        }
        
        for (int i = 0; i < xiOdomVec.size(); i+= xiOdomVec.size() - 1)
        {
            Mat8u img(480, 640);
            img.setTo(255);
            fill(img.data, img.data + 640, 0);
            fill(img.data + 640 * 479, img.data + 640 * 480, 0);
            Transf xiCamBoard = xiOdomVec[i].compose(xiCam).inverseCompose(xiBoard);
            Vector3dVec boardCam;
            xiCamBoard.transform(board, boardCam);
            Vector2dVec projectedBoard;
            cam.projectPointCloud(boardCam, projectedBoard);
            
            for (int j = 0; j < projectedBoard.size(); j++)
            {
                Vector2i p = round(projectedBoard[j]);
                img(p[1], p[0]) = 0;
            }
            
//            imshow("points", img);
            imwrite("points" + to_string(trajIdx) + to_string(i) + ".png", img);
//            waitKey();
        }         
    }
    myfile.close();
    
    
      
                  
    return 0;
    //////////////////////////////////
    //optimize the trajectory (gradient descent)
    //////////////////////////////////
    
//    Transf xiCam(0.2, 0, 0.3, 1.2, 1.2, 1.2);
    /*TrajectoryQuality * costFunction = new TrajectoryQuality(
                                            traj,
                                            xiCam,
                                            Matrix6d::Identity()*1e-3,
                                            Matrix6d::Identity());
                                            
    
    //Improvised gradient descend 
    vector<double> gradVec(paramVec.size());
    vector<double> paramVec2 = paramVec;
    double val, val2 = 999999;
    double lambda = 1e-2;
    for (int iter = 0; iter < 10000; iter++)
    {
        costFunction->Evaluate(paramVec.data(), &val, gradVec.data());
        cout << iter << " " << val <<  " " << lambda <<endl;
        double gradNorm = 0;
        for (int i = 0; i < paramVec.size(); i++)
        {
            gradNorm += gradVec[i]*gradVec[i];
        }
        if (val > val2)
        {
            lambda = max(1e-8, lambda / 3); 
            paramVec = paramVec2;
            continue;
        }
        
        val2 = val;
        paramVec2 = paramVec;
        
        if (gradNorm > lambda*lambda)
        {
            gradNorm = sqrt(gradNorm);
            for (int i = 0; i < paramVec.size(); i++)
            {
                gradVec[i] *= lambda / gradNorm;
            }
        }
        for (int i = 0; i < paramVec.size(); i++)
        {
            paramVec[i] -= gradVec[i] * lambda;
        }
    }

    for (int i = 0; i < circleCount; i++)
    {
        cout << paramVec[i * PARAM_NUM] << "   " 
            << paramVec[i * PARAM_NUM + 1] << "   " 
            << paramVec[i * PARAM_NUM + 2] << endl;
    }
    
    traj->compute(paramVec.data(), xiOdomVec, covOdomVec);
//    
//    for (int i = 0; i < xiOdomVec.size(); i++)
//    {
//        cout << xiOdomVec[i] << endl;
//        cout << covOdomVec[i].diagonal().transpose() << endl << endl;
//    }
    */
}
