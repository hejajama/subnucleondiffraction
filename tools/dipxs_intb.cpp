//
//  dipxs_intb.cpp
//  SubNucleon diffraction
//
//  Created by Heikki Mantysaari on 12/22/15.
//  Copyright © 2015 Heikki. All rights reserved.
//

#include "dipxs_intb.hpp"
#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"

struct inthelper
{
    IPGlasma *dipole;
    double r;
    double xbj;
};

double inthelperf_r(double r, void* p);
double inthelperf_theta(double theta, void* p);

int main(int argc, char* argv[])
{
    IPGlasma dipole(argv[1]);
    
    // Calculate dipole amplitude averaged over b region
    // Can be used to compare if different configurations give eventually the same t spectra
    inthelper par;
    par.dipole = &dipole;
    
    

    return 0;
}