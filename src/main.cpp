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
#include <gsl/gsl_errno.h>

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
#include "dis.hpp"
#include "virtual_photon.hpp"

using namespace std;

string InfoStr();

DipoleAmplitude *amp;

gsl_rng* global_rng;

enum MODE
{
    AMPLITUDE_DT,   // Calculate amplitude as a function of t, real/imag part set by other cli argument
    CORRECTIONS ,    // Caluclate corrections, require dipole amplitude with rotational symmetry
    PRINT_NUCLEUS,
    F2,             // Calculate F2
    SATURATION_SCALE    // Print saturation scale on a g gripd
};

int MCpoints(double t);  // Automatic mc points


int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler); // Do not let code to crash if an error occurs, we should
    // have error handling everywhere!
    
    double Qsqr=0;
    double xbj=0; // x for F2
    double t=0.1;
    int A=1;
    DGLAPDist *gd=0;  // Initialized and used if we have nucleus consisting of ipsatnucleons
    
    //double xpom=0.000959089;
    double w = 100;
    bool skewedness = false;
    double qsfluct_sigma=0;
    Fluctuation_shape fluctshape = LOCAL_FLUCTUATIONS;
    bool auto_mcintpoints = false;
    std::string wavef_file = "gauss-boosted.dat";
    
    
    cout << "# SubNucleon Diffraction by H. Mäntysaari <mantysaari@bnl.gov>, 2015-2016" << endl;
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    
    if (string(argv[1])=="-help")
    {
        cout << "-Q2, -W: set kinematics" << endl;
        cout << "-real, -imag: set real/imaginary part" << endl;
        cout << "-dipole A [ipglasma,ipsatproton] [ipglasmafile, ipsat_proton_width ipsat_proton_quark_width] [fluxtube tunbe_normalization] [albacete]" << endl;
        cout << "-corrections: calculate correction R_g^2(1+\beta^2) as a function of t. Requires rot. sym. dipole amplitude." << endl;
        cout << "-mcintpoints points/auto" << endl;
        cout << "-skewedness: enable skewedness in dipole amplitude" << endl;
        cout << "-qsfluct sigma: set width of Q_s fluctuations (0: disable); only for ipsatproton!" << endl;
        cout << "-qsfluctshape [local,quarks]: set Q_s^2 to fluctuate at each point / for each quark" << endl;
        cout << "-satscale: print saturation scale" << endl;
        cout << "-F2 Qsqr x: calculate structure function" << endl;
        cout << "-wavef_file filename" << endl;
        return 0;
    }
    
    MODE mode = AMPLITUDE_DT;
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-mcintpoints")
        {
            if (argc > i+1)
            {
                if (string(argv[i+1])=="auto")
                    auto_mcintpoints = true;
            }
            if (!auto_mcintpoints)
                MCINTPOINTS = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-Q2")
            Qsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-W")
            w=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-real")
            REAL_PART = true;
        else if (string(argv[i])=="-imag")
            REAL_PART = false;
        else if (string(argv[i])=="-dipole")
        {
            A = StrToInt(argv[i+1]);
            if (A==1)
            {
                if (string(argv[i+2])=="ipsatproton")
                {
                    amp = new Ipsat_Proton;
                    ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+3]));
                    ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+4]));
                    if (string(argv[i+5]) == "ALBACETE")
                        ((Ipsat_Proton*)amp)->SetShape(ALBACETE);
                    else
                    {
                        ((Ipsat_Proton*)amp)->SetShape(GAUSSIAN);
                        if (argc > i+5)
                        {
                            if (string(argv[i+5])=="fluxtube")
                            {
                                ((Ipsat_Proton*)amp)->SetStructure(CENTER_TUBES);
                                ((Ipsat_Proton*)amp)->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                            }
                            else if (string(argv[i+5]).substr(0,1)!="-")
                            {
                                cerr << "Unknown ipsatproton option " << argv[i+4] << endl;
                                exit(1);
                            }
                        }
                    }
                }
                else if (string(argv[i+2])=="ipglasma")
                    amp = new IPGlasma(argv[i+3]);
                else
                {
                    cerr << "Unknown dipole " << argv[i+1] << endl;
                    return -1;
                }
            }
            else
            {
                // Construct nucleus
                std::vector<DipoleAmplitude* > nucleons;
                for (int j=0; j<A; j++)
                {
                    if (string(argv[i+2])=="ipsatproton")
                    {
                        if (j==0)
                            gd = new DGLAPDist;
                        //Ipsat_Proton *nucleon = new Ipsat_Proton(gd);
                        Ipsat_Proton *nucleon = new Ipsat_Proton();
                        nucleon->SetProtonWidth(StrToReal(argv[i+3]));
                        nucleon->SetQuarkWidth(StrToReal(argv[i+4]));
                        if (string(argv[i+5]) == "ALBACETE")
                            nucleon->SetShape(ALBACETE);
                        else
                        {
                            nucleon->SetShape(GAUSSIAN);
                            if (argc > i+5)
                            {
                                if (string(argv[i+5])=="fluxtube")
                                {
                                    nucleon->SetStructure(CENTER_TUBES);
                                    nucleon->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                                }
                                else if (string(argv[i+5]).substr(0,1)!="-")
                                {
                                    cerr << "Unknown ipsatproton option " << argv[i+5] << endl;
                                    exit(1);
                                }
                            }
                        }
                        nucleons.push_back(nucleon);

                    }
                }
                amp = new Nucleons(nucleons);
            }
            
        }
        else if (string(argv[i])=="-print_nucleus")
        {
            mode = PRINT_NUCLEUS;
        }
        else if (string(argv[i])=="-skewedness")
            skewedness = true;
        else if (string(argv[i])=="-corrections")
            mode = CORRECTIONS;
        else if (string(argv[i])=="-qsfluct")
            qsfluct_sigma = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-qsfluctshape")
        {
            if (string(argv[i+1])=="local")
                fluctshape = LOCAL_FLUCTUATIONS;
            else if (string(argv[i+1])=="quarks")
                fluctshape = FLUCTUATE_QUARKS;
            else
            {
                cerr << "Unknown fluctuation type " << argv[i+1] << endl;
                exit(1);
            }
        }
        else if (string(argv[i])=="-satscale")
            mode = SATURATION_SCALE;
        else if (string(argv[i])=="-wavef_file")
            wavef_file = argv[i+1];
        else if (string(argv[i])=="-F2")
        {
            mode = F2;
            Qsqr = StrToReal(argv[i+1]);
            xbj=StrToReal(argv[i+2]);
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
        
    }
    
    
    
    
    
    // Initialize global random number generator
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);

    
    
    BoostedGauss wavef(wavef_file);
    
    
    
    
    amp->SetSkewedness(skewedness);
    if (qsfluct_sigma > 0)
    {
        ((Ipsat_Proton*)amp)->SetQsFluctuation(qsfluct_sigma);
        ((Ipsat_Proton*)amp)->SetFluctuationShape(fluctshape);
        if (A>1)
        {
            cerr << "Q_s fluctuations not implemented for nucleus!" << endl;
            exit(1);
        }
    }
    
    amp->InitializeTarget();
    

    Diffraction diff(*amp, wavef);
    
    
    cout << "# " << InfoStr() << endl;
    cout << "# " << wavef << endl;
    
    double mp = 0.938;
    double mjpsi = wavef.MesonMass();
    
    
    if (mode == PRINT_NUCLEUS)
    {
        
        
        double origin[2]={0,0};
        double max = ((IPGlasma*)amp)->MaxX();
        double min = ((IPGlasma*)amp)->MinX();
        double step =((IPGlasma*)amp)->XStep();
        cout << "# Grid min " << min << " 1/GeV, max " << max << " 1/GeV, step " << step << " 1/GeV" << endl;
        cout << "# 1/Nc(1-Tr[V(0)V(x,y)]) (re im)  1/Nc(1-Tr[V(x,y)V(x,y)])  1/Nc(Tr[1-V(x,y)])  " << endl;
        for (double y=min+step/2; y < max-step/2; y+=step)
        {
             for (double x=min+step/2; x < max-step/2; x+=step)
            {
                double p[2] = {x,y};
                
                WilsonLine &wl =((IPGlasma*)amp)->GetWilsonLine(x,y);
                double tr = wl.Trace().real();
             
                cout << y/5.068 << " " << x/5.068 << " " << ((IPGlasma*)amp)->Amplitude(0.01, origin, p) << " " << ((IPGlasma*)amp)->AmplitudeImaginaryPart(0.01, origin, p) << " " << ((IPGlasma*)amp)->Amplitude(0.01, p, p) << " " << 1.0 - tr/3.0 <<endl;
            }
         cout << endl;
        }
       
        
      /*  
        double origin[2]={0,0};
        double max = 8;
        double min = -8;
        double step = 0.1;
        cout << "# x y N(0,(x,y)) T(b) " << endl;
        for (double y=min+step/2; y < max-step/2; y+=step)
        {
            for (double x=min+step/2; x < max-step/2; x+=step)
            {
                double p[2] = {x,y};
                
                cout << y << " " << x << " " << amp->Amplitude(0.001, origin, p) << " " << ((Ipsat_Proton*)amp)->Density(Vec(x,y)) << endl;
            }
            cout << endl;
        }
        
        
        */
         
        return 0;
    }
    
    else if (mode == SATURATION_SCALE)
    {
        double max=5; int points = 100;
        for (double y=-max; y<=max; y+=2.0*max/(points-1.0))
        {
            for (double x=-max; x<=max; x+=2.0*max/(points-1.0))
            {
                //if (x*x + y*y < 0.5*0.5*5.068*5.068) // Limit maxr
                    cout << y << " " << x << "  " << amp->SaturationScale(0.001, Vec(x,y)) << endl;
            }
            cout << endl;
        }
        return 0;
    }
    
    else if (mode == AMPLITUDE_DT)
    {
        cout << "# Amplitude as a function of t, Q^2=" << Qsqr << ", W=" << w << endl;
        cout << "# t  dsigma/dt [GeV^-4] Transverse Longitudinal  " << endl;
        double tstep = 0.025;
        for (t=0; t<=2; t+=tstep)
        {
            double xpom = (mjpsi*mjpsi+Qsqr+t)/(w*w+Qsqr-mp*mp);
            if (xpom > 0.01)
            {
                cerr << "xpom = " << xpom << ", can't do this!" << endl;
                continue;
            }
            
            if(auto_mcintpoints)
                MCINTPOINTS = MCpoints(t);
            
            cout.precision(5);
            double trans = diff.ScatteringAmplitude(xpom, Qsqr, t, T);
            double lng = 0;
            if (Qsqr > 0)
                lng = diff.ScatteringAmplitude(xpom, Qsqr, t, L);

            cout << t << " ";
            cout.precision(10);
            cout << trans  << " " << lng << endl;
            
            if (t > 0.5)
                tstep = 0.1;

        }
    }
    else if (mode == CORRECTIONS)
    {
        cout << "# Real part correction" << endl;
        cout << "# t  transverse  longitudinal" << endl;
        double tstep=0.02;
        for (t=0; t<=2.5; t+=tstep)
        {
            double xpom = (mjpsi*mjpsi+Qsqr+t)/(w*w+Qsqr-mp*mp);
            if (xpom > 0.01)
            {
                cerr << "xpom = " << xpom << ", can't do this!" << endl;
                continue;
            }
            
            cout.precision(5);
            double res_t = diff.Correction(xpom, Qsqr, t, T);
            double res_l=0;
            if (Qsqr > 0)
                res_l= diff.Correction(xpom, Qsqr, t, L);
            cout << t << " ";
            cout.precision(10);
            cout << res_t << " " << res_l   << endl;
            //if (t>0.5)
            //    tstep=0.1;
        }
    }
    
    else if (mode == F2)
    {
        cout << "#F2(Qsqr=" << Qsqr << ", xbj=" << xbj << "): light charm sum " << endl;
        WaveFunction * photon = new VirtualPhoton();;
        ((VirtualPhoton*)photon)->SetQuark(Amplitude::LIGHT);
        cout << "# Quarks: " << ((VirtualPhoton*)photon)->GetParamString() << endl;
        
        amp->SetSkewedness(false);
        Diffraction f2(*amp, *photon);
        
        // Use the fact that photon-proton cross section is just diffractive amplitude at t=0
        // Note* 4pi, as convention in BoostedGaussian and VirtualPhoton classes are different!!!
        double xs_t = 4.0*M_PI*f2.ScatteringAmplitude(xbj, Qsqr, 0, T);
        double xs_l = 4.0*M_PI*f2.ScatteringAmplitude(xbj, Qsqr, 0, L);
        double structurefun = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
        
        // heavy quark contribution
        ((VirtualPhoton*)photon)->SetQuark(Amplitude::C, 1.4);
        cout << "# Quarks: " << ((VirtualPhoton*)photon)->GetParamString() << endl;
        double xs_t_c = 4.0*M_PI*f2.ScatteringAmplitude(xbj, Qsqr, 0, T);
        double xs_l_c = 4.0*M_PI*f2.ScatteringAmplitude(xbj, Qsqr, 0, L);
        double structurefun_c = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_c+xs_t_c);
        
	
        cout << structurefun << " " << structurefun_c << " " << structurefun + structurefun_c << endl;
        //DIS dis(amp);
        //cout << dis.F2(Qsqr, xbj) << endl;
        
        delete photon;
    }
    
    
    gsl_rng_free(global_rng);


    
    delete amp;
    
    if (gd != 0)
        delete gd;
    
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
    
    info << endl << amp->InfoStr();

    
    if (REAL_PART) info << "# Real part";
    else info << "# Imaginary part";
    
    if (FACTORIZE_ZINT)
        info <<". z integral factorized";
    else info << ". z integral not factorized";
    
    
    return info.str();

}

int MCpoints(double t)
{
    if (t<1)
        return 2e7;
    else if (t<2)
        return 5e7;
    else
        return 1e8;
}
