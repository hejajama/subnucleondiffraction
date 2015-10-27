/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include <iostream>
#include <vector>
#include "dipole.hpp"
#include "smooth_ws_nuke.hpp"
#include "diffraction.hpp"
#include "gauss_boost.hpp"
#include "ipsat_nucleons.hpp"
#include "vector.hpp"
using namespace std;

int main()
{
    double Qsqr=0;
    double t=0.1;
    double xpom=0.001;
    
    BoostedGauss wavef("gauss-boosted.dat");
    Smooth_ws_nuke target;
    
    Ipsat_Nucleons ipsatnuke;
    ipsatnuke.InitializeTarget();
    
    Diffraction diff(target, wavef);
    //Diffraction diff(ipsatnuke, wavef);
    
    for (t=0.01; t<=0.3; t+=0.02)
        cout << t << " " << diff.TotalCrossSection(xpom, Qsqr, t) << endl;
    
    // Try nucleus
    
    
    
    /*
    vector<Vec> nukes = ipsatnuke.GetNucleons();
    for (int i = 0; i < nukes.size(); i++)
    {
        cout << nukes[i].GetX() << " " << nukes[i].GetY() << endl;
    }
    */
    
}