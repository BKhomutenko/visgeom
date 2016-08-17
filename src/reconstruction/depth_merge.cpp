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
Depth filter algorithm. 
NOTE:
(u, v) is an image point 
(x, y) is a depth map point
*/

#include "reconstruction/depth_map.h"

#include "io.h"
#include "std.h"
#include "eigen.h"

const int COST_CHANGE = 10;
const int MAX_COST = 35; //TODO link to DEFAULT_COST

bool match(double v1, double s1, double v2, double s2)
{
    double delta = abs(v1 - v2);
    return delta < 2 * s1 or  delta < 2 * s2;
}

void filter(double & v1, double & s1, double v2, double s2)
{
    double denom = s1 + s2;
    v1 = (v1 * s2 + v2 * s1) / denom;
    s1 = s1 * s2 / denom;
}

void DepthMap::merge(const DepthMap & depth2, const Transformation<double> T12)
{
    assert((ScaleParams)(*this) == (ScaleParams)depth2);
    for (int y = 0; y < yMax; y++)
    {
        for (int x = 0; x < xMax; x++)
        {
            vector<bool> matchedMeasurementVec(depth2.hMax, false);
            //first match real hypotheses
            for (int h1 = 0; h1 < hMax; h1++)
            {
                double & cost1 = cost(x, y, h1);
                if (cost1 > MAX_COST) continue;
                bool improved = false;
                for (int h2 = 0; h2 < depth2.hMax; h2++)
                {
                    
                    if (not matchedMeasurementVec(h2) and
                            match(at(x, y, h1), sigma(x, y, h1),
                             depth2.at(x, y, h2), depth2.sigma(x, y, h2))
                    {
                        filter(at(x, y, h1), sigma(x, y, h1),
                             depth2.at(x, y, h2), depth2.sigma(x, y, h2);
                       
                        cost1 = max(cost1 - COST_CHANGE, 0);
                        improved = true;
                        matchedMeasurementVec(h2) = true;
                        break;
                    }
                }
                if (not improved) cost(x, y, h1) += COST_CHANGE;
            }
            // discard bad hypotheses
            // TODO rearrange hypotheses
            // TODO use pushHypothesis (with an overload) to add them
            int h2 = 0;
            for (int h1 = 0; h1 < hMax; h1++)
            {
                if (cost(x, y, h1) > MAX_COST)
                {
                    while (h2 < depth2.hMax and matchedMeasurementVec(h2)) h2++;
                    if (h2 < depth2.hMax)
                    {
                        at(x, y, h1) = depth2.at(x, y, h2);
                        sigma(x, y, h1) = depth2.sigma(x, y, h2);
                        //it gives the hyp 2 turns
                        cost(x, y, h1) = MAX_COST - COST_CHANGE;
                    }
                    else
                    {
                        at(x, y, h1) = DEFAULT_DEPTH;
                        sigma(x, y, h1) = DEFAULT_SIGMA;
                        cost(x, y, h1) = DEFAULT_COST;
                    }
                }
            }
        }
    }
    
    
}

