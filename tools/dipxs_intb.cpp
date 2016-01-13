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
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>

using namespace std;

gsl_rng* global_rng;
struct inthelper
{
    DipoleAmplitude *dipole;
    double r;
    double theta_r;
    double b;
    double xbj;
};

double inthelperf_b(double b, void* p);
double inthelperf_theta(double theta, void* p);

int main(int argc, char* argv[])
{
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    IPGlasma dipole(argv[1]);
    
    /*Ipsat_Proton proton;
    proton.SetProtonWidth(2);
    proton.SetQuarkWidth(2);
    proton.InitializeTarget();
    */
    
    // Calculate dipole amplitude averaged over b region
    // Can be used to compare if different configurations give eventually the same t spectra
    inthelper par;
    par.dipole = &dipole;
    par.xbj=0.01;
    par.r=1;
    par.theta_r = 0; // Dipole orientation
    
    double minb=0;
    double maxb = 2;
    
    gsl_function F;
    F.function = &inthelperf_b;
    F.params = &par;
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    
    double result,err;
    gsl_integration_qags( &F, minb, maxb, 0, 0.1, 1000, w, &result, &err );
    
    gsl_integration_workspace_free(w);
    
    double area = maxb*maxb;
    
    cout << "Average amplitude " << result/area << " +/- " << err/area << endl;
     

    return 0;
}

double inthelperf_b(double b, void* p)
{
    inthelper* par = (inthelper*)p;
    par->b = b;
    gsl_function F;
    F.function = &inthelperf_theta;
    F.params = par;
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    
    double result,err;
    gsl_integration_qags( &F, 0, 2.0*M_PI, 0, 0.1, 1000, w, &result, &err );
    
    gsl_integration_workspace_free(w);
    
    return b*result;
}

double inthelperf_theta(double theta, void* p)
{
    inthelper* par = (inthelper*)p;
    
    // quark position = b + 0.5r
    // antiquark = b - 0.5r
    
    Vec b(par->b*cos(theta), par->b*sin(theta));
    Vec rhalf(0.5*par->r*cos(par->theta_r),0.5*par->r*sin(par->theta_r) );  // 0.5*r
    
    Vec q1 = b + rhalf;
    Vec q2 = b - rhalf;
    
    double amp = par->dipole->DipoleAmplitude::Amplitude(par->xbj, q1, q2);
    double imag =par->dipole->DipoleAmplitude::AmplitudeImaginaryPart(par->xbj, q1, q2);
    
    if (std::abs(imag/amp) > 0.1)
        cerr << "Imag is important " << imag/amp << " at " << q1 << q2 << endl;
    
    return amp;
    
    
    
}
