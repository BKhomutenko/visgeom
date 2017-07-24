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

ptree serializeTransform(const Transf & xi)
{
    ptree xiNode;
    Vector3d t = xi.trans(), r = xi.rot();
    ptree valNode;
    valNode.put_value(t[0]); xiNode.push_back( make_pair( "", valNode) );
    valNode.put_value(t[1]); xiNode.push_back( make_pair( "", valNode) );
    valNode.put_value(t[2]); xiNode.push_back( make_pair( "", valNode) );
    valNode.put_value(r[0]); xiNode.push_back( make_pair( "", valNode) );
    valNode.put_value(r[1]); xiNode.push_back( make_pair( "", valNode) );
    valNode.put_value(r[2]); xiNode.push_back( make_pair( "", valNode) );
    return xiNode;
}

int main(int argc, char** argv) 
{
    
    ptree odomRoot, imgRoot;
    read_json(argv[1], odomRoot);
    read_json(argv[2], imgRoot);
    
    map<double, Transf> transfMap;
    for (auto & x : odomRoot)
    {
        double t = x.second.get<int>("time_s") + x.second.get<int>("time_ns") * 1e-9;
        Transf xi = readTransform(x.second.get_child("pose"));
        transfMap[t] = xi;
    }

    map<double, string> nameMap;
    for (auto & x : imgRoot)
    {
        double t = x.second.get<int>("time_s") + x.second.get<int>("time_ns") * 1e-9;
        string name = x.second.get<string>("fname");
        nameMap[t] = name;
    }
    
    
    ptree odomTree, imgTree;
    for (auto & x : nameMap)
    {
        auto it = transfMap.lower_bound(x.first);
        if (it == transfMap.begin())
        {
            ptree xiNode = serializeTransform(it->second);
            odomTree.push_back( make_pair("", xiNode) );
            ptree valNode;
            valNode.put_value(x.second);
            imgTree.push_back( make_pair("", valNode) );   
        }
        else if (it != transfMap.end())
        {
            auto it0 = it;
            it0--;
            double t0 = it0->first;
            double t1 = it->first;
            Transf xi0 = it0->second;
            Transf xi1 = it->second;
            Transf zeta = xi0.inverseCompose(xi1);
            double lambda = (x.first - t0) / (t1 - t0);
            zeta.scale(lambda);
            Transf xi = xi0.compose(zeta);
            ptree xiNode = serializeTransform(xi);
            odomTree.push_back( make_pair("", xiNode) );
            ptree valNode;
            valNode.put_value(x.second);
            imgTree.push_back( make_pair("", valNode) );         
        }
    }
    ptree pout;
    pout.add_child("names", imgTree);
    pout.add_child("odom", odomTree);
    write_json("prepared.json", pout);
}



