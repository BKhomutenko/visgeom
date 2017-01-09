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

class CircularTrajectory : public ITrajectory
{
public:
    CircularTrajectory(int numberSteps, double period, int numberCircles = 1) : 
        _numberSteps(numberSteps),
        _period(period),
        _numberCircles(numberCircles) {}
        
    virtual vector<Transf> compute(const double * params) const
    {
        vector<Transf> res;
        double tMax = _numberSteps * _period;
        for (int circle = 0; circle < _numberCircles; circle++)
        {
            double omega = atan(params[circle]) / tMax;
            double v = 1;
            double alpha = (omega * _period);
            double dist = v * _period;
            Transf dxi(-dist * sin(alpha*0.5), dist * cos(alpha*0.5), 0, 0, 0, alpha);
            
            res.push_back(dxi);
            for (int i = 1; i < _numberSteps; i++)
            {
                res.push_back(res.back().compose(dxi));
            }
        }
        return res;
    }
    
    virtual int paramSize() const { return _numberCircles; }
    
    int _numberSteps;
    double _period;
    int _numberCircles;
};


int main(int argc, char** argv) 
{
    int circleCount = 2;
    vector<double> paramVec;
    for (int i = 0; i < circleCount; i++)
    {
        paramVec.push_back(i * 0.1);
    }
    Transf xiCam(0.2, 0, 0.3, 1.2, 1.2, 1.2);
    double period = 0.1;
    double numberSteps = 100;
    TrajectoryQuality * costFunction = new TrajectoryQuality(
                                            new CircularTrajectory(numberSteps, period, circleCount),
                                            xiCam,
                                            Matrix6d::Identity(),
                                            Matrix6d::Identity());
                                            
    ceres::GradientProblem problem(costFunction);

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, paramVec.data(), &summary);

    cout << summary.FullReport() << "\n";
    double tMax = numberSteps * period;
    for (int i = 0; i < circleCount; i++)
    {
        cout << atan(paramVec[i]) / tMax << endl;
    }
}

