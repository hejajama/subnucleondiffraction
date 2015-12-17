/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include <gsl/gsl_rng.h>

#include <tools/tools.hpp>


#include "dipole.hpp"
#include "smooth_ws_nuke.hpp"
#include "diffraction.hpp"
#include "gauss_boost.hpp"
#include "ipsat_nucleons.hpp"
#include "ipsat_proton.hpp"
#include "vector.hpp"
#include "subnucleon_config.hpp"
#include "ipglasma.hpp"
#include "nucleons.hpp"

using namespace std;

string InfoStr();

DipoleAmplitude *amp;

gsl_rng* global_rng;

int main(int argc, char* argv[])
{
    double Qsqr=0;
    double t=0.1;
    double xpom=0.001;
    PROCESS p = COHERENT;
    bool print_nucleus = false;
    
    cout << "# SubNucleon Diffraction by H. Mäntysaari <mantysaari@bnl.gov>, 2015" << endl;
    
    if (string(argv[1])=="-help")
    {
        cout << "-real, -imag: set real/imaginary part" << endl;
        cout << "-dipole [ipsat,ipnonsat,ipglasma,ipsatproton,nucleons] [ipglasmafile, ipsat_radius_fluctuation_fraction, ipsat_proton_width ipsat_proton_quark_width]" << endl;
        cout << "-mcintpoints points" << endl;
        return 0;
    }
        
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-coherent")
            p = COHERENT;
        else if (string(argv[i])=="-incoherent")
            p = INCOHERENT;
        else if (string(argv[i])=="-mcintpoints")
            MCINTPOINTS = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-real")
            REAL_PART = true;
        else if (string(argv[i])=="-imag")
            REAL_PART = false;
        else if (string(argv[i])=="-dipole")
        {
            if (string(argv[i+1])=="ipsat" or string(argv[i+1])=="ipnonsat")
            {
                amp = new Ipsat_Nucleons;
                if (string(argv[i+1])=="ipsat")
                    ((Ipsat_Nucleons*)amp)->SetSaturation(true);
                else
                    ((Ipsat_Nucleons*)amp)->SetSaturation(false);
                ((Ipsat_Nucleons*)amp)->SetFluctuatingNucleonSize(StrToReal(argv[i+2]));
                
            }
            else if (string(argv[i+1])=="ipsatproton")
            {
                amp = new Ipsat_Proton;
                ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+2]));
                ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+3]));
                ((Ipsat_Proton*)amp)->SetShape(EXPONENTIAL);
                
            }
            else if (string(argv[i+1])=="ipglasma")
                amp = new IPGlasma(argv[i+2]);
            else if (string(argv[i+1])=="nucleons")
            {
                amp = new Nucleons;
                ((Nucleons*)amp)->SetProtonWidth(StrToReal(argv[i+2]));
                ((Nucleons*)amp)->SetQuarkWidth(StrToReal(argv[i+3]));
            }
            else
            {
                cerr << "Unknown dipole " << argv[i+1] << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-print_nucleus")
        {
            print_nucleus = true;
        }
    }
    
    
    
    
    
    //
    //IPGlasma glasma("data/V.dat");
    //IPGlasma glasma("proton_samples/proton_tarkka");
    double origin[2]={0,0};
    
    double max = ((IPGlasma*)amp)->MaxX();
    double min = ((IPGlasma*)amp)->MinX();
    double step =((IPGlasma*)amp)->XStep();

    for (double y=min+step/2; y < max-step/2; y+=step)
    //for (double y=-4.8; y < 4.8; y+=0.2)
    {
        for (double x=min+step/2; x < max-step/2; x+=step)
        //for (double x=-4.8; x < 4.8; x+=0.2)
        {
            //origin[0]=x; origin[1]=y;
            double p[2] = {x,y};
            
            cout << y << " " << x << " " << ((IPGlasma*)amp)->Amplitude(0.01, origin, p) << " " << ((IPGlasma*)amp)->Amplitude(0.01, p, p) << endl;
        }
        cout << endl;
       
        
    }
     
     return 0;
   
    
    // Initialize global random number generator
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);

    
    
    BoostedGauss wavef("gauss-boosted.dat");
    
    amp->InitializeTarget();
    

    Diffraction diff(*amp, wavef);
    
    
    cout << "# " << InfoStr() << endl;
    cout << "# " << wavef << endl;
    
    
    if (print_nucleus)
    {
        // Print ipsat nucleus, todo: ipglasma
        
            std::vector<Vec> positions = ((Ipsat_Proton*)amp)->GetQuarks();
            for (int j=0; j<3; j++)
                cout << positions[j].Len() << endl;
        
        return 0;
        //std::vector<Vec> positions = ((Ipsat_Proton*)amp)->GetQuarks();
        std::vector<double> radii =((Ipsat_Proton*)amp)->GetRadii();
        cout << "# x   y    radius   [GeV^-1]" << endl;
        for (int i=0; i<positions.size(); i++)
        {
            cout << positions[i].GetX() << " " << positions[i].GetY() << " " << radii[i] << endl;
        }

        /*
        std::vector<Vec> positions = ((Ipsat_Nucleons*)amp)->GetNucleons();
        std::vector<double> Bps = ((Ipsat_Nucleons*)amp)->GetB_ps() ;
        cout << "# x   y    radius   [GeV^-1]" << endl;
        for (int i=0; i<positions.size(); i++)
        {
            cout << positions[i].GetX() << " " << positions[i].GetY() << " " << std::sqrt(2.0*Bps[i]) << endl;
        }*/
        return 0;
        
    }
    
    if (p == INCOHERENT)
        cout << "# t    dsigma/dt [GeV^4] " << endl;
    if (p == COHERENT)
        cout << "# t    Re or Im A [GeV^2] " << endl;
    for (t=0.0; t<=2.61; t+=0.150)
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
    
    
    gsl_rng_free(global_rng);


    
    delete amp;
    
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
    
    info << " Dipole: " << amp->InfoStr();

    
    if (REAL_PART) info << ". Real part";
    else info << ". Imaginary part";
    
    if (FACTORIZE_ZINT)
        info <<". z integral factorized";
    else info << ". z integral not factorized";
    
    info << ". Corrections: ";
    if (CORRECTIONS) info << "enabled";
    else info << "disabled";
    
    return info.str();

}