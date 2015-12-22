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
    //Ipsat_Proton dipole;
    //dipole.SetProtonWidth(0.0001);
    //dipole.SetQuarkWidth(4);
    //dipole.InitializeTarget();
    //cout << "initialized\n";
    
    for (double y=-6; y<6; y+=0.05)
    {
        for (double x=-6; x<6; x+=0.05)
        {
            cout << dipole.SaturationScale(0.01, Vec(0,0)) << " ";
        }
        cout << endl;
    }
    
    
    gsl_rng_free(global_rng);
    return 0;
}  