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

#include "std.h"
#include "json.h"
#include "io.h"
#include "geometry/geometry.h"

int main(int argc, char** argv) 
{
    
    ptree root1, root2;
    read_json(argv[1], root1);
    read_json(argv[2], root2);
    
    map<double, string> nameMap1, nameMap2;
    for (auto & x : root1)
    {
        double t = x.second.get<int>("time_s") + x.second.get<int>("time_ns") * 1e-9;
        string name = x.second.get<string>("fname");
        nameMap1[t] = name;
    }
    for (auto & x : root2)
    {
        double t = x.second.get<int>("time_s") + x.second.get<int>("time_ns") * 1e-9;
        string name = x.second.get<string>("fname");
        nameMap2[t] = name;
    }
    
    const double SYNC_PRESICION = 0.01;
    ptree tree1, tree2;
    for (auto & x : nameMap1)
    {
        auto it = nameMap2.lower_bound(x.first);
        if (it != nameMap2.begin() and it != nameMap2.end())
        {
            auto it0 = it;
            it0--;
            double t0 = it0->first;
            double t1 = it->first;
            double t = x.first;
            if (min(abs(t - t1), abs(t - t0)) < SYNC_PRESICION)
            {
                ptree valNode;
                valNode.put_value(x.second);
                tree1.push_back( make_pair("", valNode));  
                if (abs(t - t1) < abs(t - t0)) // it points to the matching image
                {
                    valNode.put_value(it->second);
                    tree2.push_back( make_pair("", valNode));  
                }
                else
                {
                    valNode.put_value(it0->second);
                    tree2.push_back( make_pair("", valNode));  
                }
            }
        }
    }
    ptree pout;
    pout.add_child("img1", tree1);
    pout.add_child("img2", tree2);
    write_json("prepared.json", pout);
}



