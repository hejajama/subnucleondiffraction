//
//  satscale.cpp
//  SubNucleon diffraction
//  2015 Heikki Mäntysaari
//

#include <iostream>
#include <gsl/gsl_rng.h>
#include "../src/ipglasma.hpp"
#include "../src/vector.hpp"
#include "../src/ipsat_proton.hpp"

using namespace std;
gsl_rng* global_rng;
int main(int argc, char* argv[])
{
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    IPGlasma dipole(argv[1]);
    
    Ipsat_Proton dipole_ipsat;
    dipole_ipsat.SetProtonWidth(0.0001);
    dipole_ipsat.SetQuarkWidth(4);
    dipole_ipsat.InitializeTarget();
   
    
    //cout << "initialized\n";
    
    double avgscale=0;
    int points=0;
    
    for (double y=-6; y<6; y+=0.05)
    {
        for (double x=-6; x<6; x+=0.05)
        {
            double qs = dipole.SaturationScale(0.01, Vec(x,y));
            if (qs < 0)
            {
                qs=0;
                //cerr << "Cant solve satscale at point " << x << ", " << y << endl;
                continue;
            }
            avgscale += qs; points++;
            //cout << y << " " << x << " " << qs << endl;
        }
        //cout << endl;
    }
    if (points < 100)
        cerr << "NOTE: only " << points << " obtained" << endl;
    cout << "#AVG satscale " << avgscale / points << endl;
    
    
    gsl_rng_free(global_rng);
    return 0;
}  