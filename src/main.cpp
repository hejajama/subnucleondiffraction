/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "dipole.hpp"
#include "smooth_ws_nuke.hpp"
#include "diffraction.hpp"
#include "gauss_boost.hpp"
#include "ipsat_nucleons.hpp"
#include "vector.hpp"
#include "subnucleon_config.hpp"
#include "ipglasma.hpp"

using namespace std;

string InfoStr();

bool saturation = true; // For ipsat


int main(int argc, char* argv[])
{
    double Qsqr=0;
    double t=0.1;
    double xpom=0.001;
    PROCESS p = COHERENT;
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-coherent")
            p = COHERENT;
        else if (string(argv[i])=="-incoherent")
            p = INCOHERENT;
        else if (string(argv[i])=="-nonsat")
            saturation = false;
        
    }
    
    
    //
    //IPGlasma glasma("data/V.dat");
    
    //double origin[2]={0,0};
    
    /*
    for (double y=-11.5; y < 11.5; y+=0.05)
    {
        for (double x=-11.5; x < 11.5; x+=0.05)
        {
            double p[2] = {x,y};
            cout << y << " " << x << " " << glasma.Amplitude(0.01, origin, p) << endl;
        }
        cout << endl;
        
    }*/
   
    
    
    BoostedGauss wavef("gauss-boosted.dat");
    Smooth_ws_nuke target;
    
    Ipsat_Nucleons ipsatnuke;
    ipsatnuke.InitializeTarget();
    ipsatnuke.SetSaturation(saturation);
    
    //IPGlasma glasma("data/V.dat");
    
    //Diffraction diff(target, wavef);
    Diffraction diff(ipsatnuke, wavef);
    //Diffraction diff(glasma, wavef);
    
    cout << "# SubNucleon Diffraction" << endl;
    cout << "# " << InfoStr() << endl;
    cout << "# " << wavef << endl;
    
    if (p == INCOHERENT)
        cout << "# t    dsigma/dt [GeV^4] " << endl;
    if (p == COHERENT)
        cout << "# t    Re A [GeV^2] " << endl;
    for (t=0.0; t<=0.21; t+=0.005)
    {
        double res = 0;
        cout.precision(5);
        if (p == INCOHERENT)
            res =diff.TotalCrossSection(xpom, Qsqr, t);
        else if (p == COHERENT)
            res = diff.CoherentCrossSection(xpom, Qsqr, t);
        cout << t << " ";
        cout.precision(10);
        cout << res  << endl;

    }
    
    // Try nucleus
    
    
    
    /*
    vector<Vec> nukes = ipsatnuke.GetNucleons();
    for (int i = 0; i < nukes.size(); i++)
    {
        cout << nukes[i].GetX() << " " << nukes[i].GetY() << endl;
    }
    */
    
}


string InfoStr()
{
    stringstream info;
    
    info << "Parameters: MCINTPOINTS: " << MCINTPOINTS << " ZINT_INTERVALS " << ZINT_INTERVALS << " MCINTACCURACY " << MCINTACCURACY << " ZINT _RELACCURACY " << ZINT_RELACCURACY;
    info << ". Integration method ";
    if (MCINT == MISER)
        info << "MISER";
    else if (MCINT == VEGAS)
        info << "VEGAS";
    else
        info << "unknown!";
    info <<" Saturation: ";
    if (saturation)
        info << "enabled";
    else
        info << "disabled";
    
    return info.str();

}