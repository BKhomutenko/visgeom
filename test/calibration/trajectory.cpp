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

const int PARAM_NUM = 4;
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
        
        double dist = params[2];
        double alpha = params[3];
        double ca = cos(alpha*0.5);
        double sa = sin(alpha*0.5);
        Transf xi0(0, params[0], 0,
                     0, 0, params[1]);
        //TODO modelizer proprement l'odometrie            
        covVW <<    1e-4,        0, 
                    0,       1e-4;
                
        //circular motion model
        Transf dxi(-dist * sa, dist * ca, 0, 0, 0, alpha);
        //motion jacobian
        Matrixd<6, 2> dxidu;
        dxidu <<    -sa,     -dist/2 * ca,
                    ca,     -dist/2 * sa,
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

void fourParamsAnalysis(const double * const params, const Transf xiBoard)
{
    cout << "Parameter analysis" << endl;
    
    const double R = params[3] / params[0];
    const double D = -params[2] + xiBoard.trans()[0];
    const double th = params[1];
    cout << "R : " << R << endl;
    cout << "D : " << D << endl;
    cout << "th : " << th << " or  " << th / M_PI * 180 << " deg" << endl;
    
    
    const double x = D + R * sin(th);
    const double y = R * cos(th);
    
    cout << "center coordinates : " <<  x << "   " << y << endl;
    
}

int main(int argc, char** argv) 
{

//    vector<Transf> subVec {Transf( 3.50115, 0.349983, 0.499947, -1.21177,  1.20848, -1.20952),
//                                Transf( 3.56777 ,0.360773 ,0.490473 , -1.24337 , 1.19484, -1.24086),
//                                Transf( 3.76132,  0.41007 ,0.464887 , -1.32837,  1.15668, -1.32002),
//                                Transf( 4.03562 , 0.49536 ,0.432115 , -1.47058 ,  1.0778 ,-1.45874),
//                                Transf( 4.3787, 0.614852, 0.402354 , -1.65827, 0.941195, -1.65902)};

     /*   vector<Transf> subVec {Transf( 3.48345, 0.34885  , 0 , -1.20864,  1.20987,  -1.21105),
                                Transf( 3.53387,  0.32341 , 0 , -1.21527,  1.21199,  -1.2225 ),
                                Transf( 3.67889, 0.262883 , 0 , -1.23466,  1.21828,  -1.25445),
                                Transf( 3.89535, 0.166555 , 0 , -1.26156,  1.22972,  -1.30231),
                                Transf( 4.16502, 0.0273281,  0, -1.28958,  1.24757,  -1.36067)};
    vector<Transf> optimVec {Transf( 3.49097, 0.34758, 0,  -1.20894,  1.20852, -1.20934),
                                Transf( 3.49044, 0.346785, 0,  -1.21047,  1.20729, -1.20665),
                                Transf( 3.51464, 0.342563, 0,   -1.2163,  1.20731, -1.20437),
                                Transf( 3.56106, 0.335435, 0,  -1.22571,  1.20809, -1.20326),
                                Transf( 3.62818, 0.328741,0,  -1.23619,   1.2091, -1.20371)};
                                
 

//vector<Transf> stereoVec {Transf( 3.48284, 0.348637 ,0 , -1.21019  ,1.20978,  -1.2083),
//                                Transf( 3.49206,  0.356478, 0,  -1.21078 , 1.21011, -1.20566),
//                                Transf( 3.51522,  0.358268, 0,  -1.21328 , 1.21054, -1.20487),
//                                Transf( 3.55623,  0.355995, 0,  -1.21787,  1.21059, -1.20602),
//                                Transf( 3.61632, 0.352219 ,0 , -1.22372,  1.21071, -1.20881)};
    Matrix3d RGT;
    RGT <<    0,      0,      1, 
            -1,     0,      0,
            0,      -1,     0;
    Transf xiGT(Vector3d(3.5, 0.35, 0), RGT);
    cout << xiGT << endl;
    for (int i = 0; i < 5; i++)
    {
        Transf err1 = xiGT.inverseCompose(optimVec[i]);
        Transf err2 = xiGT.inverseCompose(subVec[i]);
        cout << i * 0.002 << " & " << i * 0.01  << std::setprecision(3)
                << " & " << err1.trans().norm() << " & " << err1.rot().norm()
                << " & " << err2.trans().norm() << " & " << err2.rot().norm()
                << "\\\\ \\hline" << endl; 
    }
    
    return 0;*/
    
    ptree root;
    read_json(argv[1], root);
    
    
    int circleCount = 2;
    vector<ITrajectory*> trajVec;
    double numberSteps = 30;
    vector<double> paramVec;
    
    paramVec.push_back(root.get<double>("trajectory1.y_0"));
    paramVec.push_back(root.get<double>("trajectory1.theta_0"));
    paramVec.push_back(root.get<double>("trajectory1.v"));
    paramVec.push_back(root.get<double>("trajectory1.w"));
    trajVec.push_back(new CircularTrajectory(root.get<int>("trajectory1.number_steps")));

    paramVec.push_back(root.get<double>("trajectory2.y_0"));
    paramVec.push_back(root.get<double>("trajectory2.theta_0"));
    paramVec.push_back(root.get<double>("trajectory2.v"));
    paramVec.push_back(root.get<double>("trajectory2.w"));
    trajVec.push_back(new CircularTrajectory(root.get<int>("trajectory2.number_steps")));
    
    vector<Transf> xiOdomVec;
    vector<Matrix6d> covOdomVec;
    
//    vector<Transf> xiOdomBiasVec;
//    vector<Matrix6d> covOdomBiasVec;
    
    
    vector<double> camParams = readVector<double>(root.get_child("camera_params"));
    const int WIDTH = root.get<int>("camera_width");
    const int HEIGHT = root.get<int>("camera_height");

    EnhancedCamera cam(WIDTH, HEIGHT, camParams.data());



    //ZOE
//    const int WIDTH = 800;
//    const int HEIGHT = 600;
//    Transf xiBoard(7, M * step / 2, N * step + 0.2, 0, 0, 0);
//    Matrix3d R;
//    R <<    0,      0,      1, 
//            -1,     0,      0,
//            0,      -1,     0;
//    Transf xiCam(Vector3d(3.5, 0.35, 0.3), R);
//    cout << xiCam.rot() << endl;
//    vector<double> camParams{0.55, 1.1, 200, 200, WIDTH / 2, HEIGHT / 2};
//    EnhancedCamera cam(WIDTH, HEIGHT, camParams.data());
//    const double curvatureMax = 1. / (2.588 / tan(M_PI / 6) + 1.);
//    
    //PIONEER
//    Transf xiBoard(3, -0.4, 0.5, 0, 0, 0);
//    Matrix3d R;
//    R <<    0,      0,      1, 
//            -1,     0,      0,
//            0,      -1,     0;
//    Transf xiCam(Vector3d(0.3, 0, 0.3), R);
//    cout << xiCam.rot() << endl;
//    vector<double> camParams{0.719981,  1.03894,  381.974,  382.378,  523.6,  366.235};
//    EnhancedCamera cam(1024, 768, camParams.data());
    
//    Transf xiCam2(0, 0, 0, 0, 0, 0);
    
    
    
    TrajectoryVisualQuality * quality = new TrajectoryVisualQuality(
                                                trajVec, root.get_child("quality_parameters"), &cam);
    
    //////////////////////////////////
    //Test the visual covariance
    //////////////////////////////////
    Transf xiCam = readTransform(root.get_child("quality_parameters.xiBaseCam"));
    Transf xiBoard = readTransform(root.get_child("quality_parameters.xiOrigBoard"));
    Matrix6d C = quality->visualCov(Transf(0.1, 0.1, 0, 0, 0, .1).compose(xiCam));
    cout << C << endl << endl;
    cout << C.inverse() << endl;
    cout << quality->EvaluateCost(paramVec.data()) << endl;
    
    
    
    //////////////////////////////////
    //Solve the optimization problem
    //////////////////////////////////

    ceres::GradientProblem problem(quality);
    ceres::GradientProblemSolver::Options options;
//    options.max_num_iterations = 25;
//    options.max_num_line_search_step_size_iterations = 100;
    options.max_num_iterations = 1000;
    options.minimizer_progress_to_stdout = true;
    ceres::GradientProblemSolver::Summary summary;
    if (root.get<bool>("optimize"))
    {
        ceres::Solve(options, problem, paramVec.data(), &summary);
        std::cout << summary.FullReport() << endl;
    }
    
    for (int i = 0; i < circleCount; i++)
    {
        for (int j = 0; j < PARAM_NUM; j++)
        {
             cout << paramVec[i * PARAM_NUM + j] << "   ";
        }
        
        cout <<  endl;
    }
    
    
    //output the trajectory parameters
    //TODO change back everything
//    fourParamsAnalysis(paramVec.data(), xiBoard);
//    fourParamsAnalysis(paramVec.data() + PARAM_NUM, xiBoard);
    
    ofstream myfile;
    myfile.open("/home/bodyk/projects/python/trajectory/traj3");
    ofstream jsonfile;
    jsonfile.open("traj3.json");
    
    ofstream jsonBiasfile;
    jsonBiasfile.open("traj3Bias.json");
    
    vector<double> paramBiasVec = paramVec;
//    paramBiasVec[0] *= 1.03;
//    paramBiasVec[4] += 0.0015;
//    paramBiasVec[3] *= 1.03;
//    paramBiasVec[7] += 0.0015;
    for (int trajIdx = 0; trajIdx < trajVec.size(); trajIdx++)
    {
        auto traj = trajVec[trajIdx];
        
        traj->compute(paramVec.data() + trajIdx * traj->paramSize(), xiOdomVec, covOdomVec);
//        traj->compute(paramBiasVec.data() + trajIdx * traj->paramSize(), xiOdomBiasVec, covOdomBiasVec);
        
        for (int i = 0; i < xiOdomVec.size(); i++)
        {
            //to read as json
            
            jsonfile << "            [" 
                        << xiOdomVec[i].trans()[0] << ", " 
                        << xiOdomVec[i].trans()[1] << ", " 
                        << xiOdomVec[i].rot()[2] << "],"<< std::endl;
//            jsonBiasfile << "            [" 
//                        << xiOdomBiasVec[i].trans()[0] << ", " 
//                        << xiOdomBiasVec[i].trans()[1] << ", " 
//                        << xiOdomBiasVec[i].rot()[2] << "],"<< std::endl;
            //to plot
            myfile << xiOdomVec[i].trans()[0] << "   " << xiOdomVec[i].trans()[1]  << endl;
            
//            cout << quality->imageLimitsCost(xiOdomVec[i].compose(xiCam)) << endl;
        }
        
        for (int i = 0; i < xiOdomVec.size(); i+= /*xiOdomVec.size() - */1)
        {
            Mat8u img(HEIGHT, WIDTH);
            img.setTo(255);
            fill(img.data, img.data + WIDTH * 2, 0);
            fill(img.data + WIDTH * (HEIGHT - 2), img.data + WIDTH * HEIGHT, 0);
            Transf xiCamBoard = xiOdomVec[i].compose(xiCam).inverseCompose(xiBoard);
            Vector3dVec boardCam;
            xiCamBoard.transform(quality->_board, boardCam);
            Vector2dVec projectedBoard;
            cam.projectPointCloud(boardCam, projectedBoard);
            
            cout << xiOdomVec[i] << "      " 
                << quality->imageLimitsCost(xiOdomVec[i].compose(xiCam)) << endl;
            for (int j = 0; j < projectedBoard.size(); j++)
            {
                Vector2i p = round(projectedBoard[j]);
                img(p[1], p[0]) = 0;
                img(p[1]-1, p[0]) = 0;
                img(p[1]+1, p[0]) = 0;
                img(p[1], p[0]-1) = 0;
                img(p[1], p[0]+1) = 0;
            }
            
            imshow("points", img);
//            imwrite("points" + to_string(trajIdx) + to_string(i) + ".png", img);
            waitKey();
        }      
    }
    myfile.close();
    jsonfile.close();
    jsonBiasfile.close();
    
      
                  
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

