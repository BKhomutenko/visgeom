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
#include "ocv.h"
#include "io.h"
#include "eigen.h"
#include "ceres.h"

void calcHist(const Mat8u & img, int valMin, int valMax,
        int binSize, vector<int> & hist, double correction = 1)
{
    int binCount = (valMax - valMin) / binSize;
    hist.resize(binCount);
    fill(hist.begin(), hist.end(), 0);
    double surplus = 0;
    for (int v = 0; v < img.rows; v++)
    {
        for (int u = 0; u < img.cols; u++)
        {   
            int val = img(v, u);
            if (val >= valMax) continue;
            double scaled = val * correction;
            int rounded = round(scaled);
            surplus = scaled - rounded;
            int idx = (rounded - valMin) / binSize;
            if (idx > 0 and idx < binCount) hist[idx]++;
        }
    }
}

int main (int argc, char const* argv[])
{
    Mat8u img1 = imread("/home/bogdan/projects/data/stereo/equalization/left-right/left01.png", 0);
    Mat8u img2 = imread("/home/bogdan/projects/data/stereo/equalization/left-right/right01.png", 0);
//    img1 = img1.rowRange(74, 363).colRange(144, 496).clone();
//    img2 = img2.rowRange(74, 363).colRange(144, 496).clone();
    GaussianBlur(img1.rowRange(74, 363).colRange(144, 496), img1, Size(5, 5), 0, 0);
    GaussianBlur(img2.rowRange(74, 363).colRange(144, 496), img2, Size(5, 5), 0, 0);
    
    resize(img1, img1, Size(0, 0), 0.5, 0.5);
    resize(img2, img2, Size(0, 0), 0.5, 0.5);
    double K = 1;
    
    vector<int> hist1, hist2;
    
    int binSize = 8;
    
    calcHist(img1, 0, 248, binSize, hist1);
    /*int errAcc1 = 0;
    for (int i = 0; i < hist1.size(); i++)
    {
        if (hist2[i] == 0) break;
        int err = hist2[i] - hist1[i];
        cout << setw(10) << hist1[i] << setw(10) << err  << setw(10) << hist2[i] << endl;
        errAcc1 += abs(err);
    }*/
/*
    // let image be scaled by alpha ~= 1
    // gamma = (alpha - 1) / alpha
    // alpha = 1 / (1 - gamma)
    vector<double> gammaVec;
    for (int i = 0; i < hist1.size() and hist2[i] != 0; i++)
    {
        double num = hist1[i] - hist2[i];
        double denom;
        if (i == 0) //edge case 1
        {
            denom = binSize * (hist2[1] - hist2[0]);
        }
        else if (i == hist1.size() - 1)  //edge case 2
        {
            // denom = binSize * (i * hist2[i - 1] - (2 * i + 1) * hist2[i]);
            continue;
        }
        else  //general case
        {
            denom = binSize * (i * hist2[i - 1]  + (i + 1) * hist2[i + 1] - (2 * i + 1) * hist2[i]);
        }
        //to avoid infty
        if (abs(denom) < 1) continue;
        cout << "num / den " << setw(10) << num << setw(10) << denom << setw(15) << num / denom << endl;
        gammaVec.push_back(num / denom);
    }
    sort(gammaVec.begin(), gammaVec.end());
    double gamma = gammaVec[gammaVec.size() / 2];
    double alpha = 1. / (1. - gamma);
    cout << "gamma / alpha : " << gamma<< setw(10) << alpha << endl;*/
    
    double minAlpha = 100;
    double minErr = 1e10;
    
    for (double alpha = 0.5; alpha < 2; alpha += 0.001)
    {
        cout << alpha << endl;   
        calcHist(img2, 0, 248, binSize, hist2, alpha);
        int errAcc2 = 0;
        for (int i = 0; i < hist1.size(); i++)
        {
            int err = hist2[i] - hist1[i];
//            cout << setw(10) << hist1[i] << setw(10) << err  << setw(10) << hist2[i] << endl;
            errAcc2 += abs(err);
        }
        if (minErr > errAcc2)
        {
            minErr = errAcc2;
            minAlpha = alpha;
        }
        cout << "error accumulator : " << setw(10) << errAcc2 << endl;
    }
    cout << "The best alpha : " << minAlpha << endl;
    return 0;
}

