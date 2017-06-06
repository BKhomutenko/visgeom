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

#include "render/texture.h"
#include "render/aa_filter_lt.h"

Texture::Texture(const Mat8u & textureMat) :
    _textureMat(textureMat.clone()),
    _grid(_textureMat.cols, _textureMat.rows, _textureMat.data),
    _interpolator(_grid) {}
    
Texture::~Texture() {}

uchar Texture::sample(const Vector2d & pt, const Matrix2d & basis) const
{
    double res = 0;
    const Vector2d eu = basis.col(0);
    const Vector2d ev = basis.col(1);
    //fixme implement a pyramid
    const LookupFilter &  filter = LookupFilter::instance();
    int nu = min(max(round(eu.norm() * filter.getRadius()), 1.), 24.); //FIXME make mipmap
    int nv = min(max(round(ev.norm() * filter.getRadius()), 1.), 24.);
    if (nu == 1 and nv == 1 or 1) 
    {
        _interpolator.Evaluate(pt[1], pt[0], &res, NULL, NULL);
//        res = bilinear<double>(_textureMat, pt[0], pt[1]);
    }
    else
    {
//        nu = max(nu, 2);
//        nv = max(nv, 2);
        Vector2d ustep = eu * filter.computeStep(nu);
        Vector2d vstep = ev * filter.computeStep(nv);
        Vector2d pt0 = pt + eu * filter.computeOrigin(nu) + ev * filter.computeOrigin(nv);
        
        //filter provides normalized kernels
        const vector<double> & uKernel = filter.getKernel(nu);
        const vector<double> & vKernel = filter.getKernel(nv);
        for (int vidx = 0; vidx < nv; vidx++)
        {
            Vector2d ptCur = pt0 + vidx * vstep;
            double ures = 0;
            for (int uidx = 0; uidx < nu; uidx++, ptCur += ustep)
            {
                double sample = bilinear<double>(_textureMat, ptCur[0], ptCur[1]);

//                double sample = 1;
//                _interpolator.Evaluate(ptCur[1], ptCur[0], &sample, NULL, NULL);
                ures += sample * uKernel[uidx] ;
            }
            res += ures * vKernel[vidx];
        }
//        res /= nu * nv;
    }
    return uchar( min(255., max(res, 0.)) );
}

