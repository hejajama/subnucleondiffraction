/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2021
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include <gsl/gsl_rng.h>

#include <gsl/gsl_errno.h>

#include "dipole.hpp"
#include "diffraction.hpp"
#include "gauss_boost.hpp"
#include "gaus_lc.h"
#include "nrqcd_wf.hpp"
#include "smooth_ws_nuke.hpp"
#include "ipsat_proton.hpp"
#include "vector.hpp"
#include "subnucleon_config.hpp"
#include "ipglasma.hpp"
#include "nucleons.hpp"
#include "dis.hpp"
#include "virtual_photon.hpp"
#include "dvcs_photon.hpp"
#include "gitsha1.h"

using namespace std;

DipoleAmplitude *amp;

gsl_rng* global_rng;

enum MODE
{
    AMPLITUDE_DT,   // Calculate amplitude as a function of t, real/imag part set by other cli argument
    TOTALCROSSSECTION,  // calculate quantities for computing total cross section
    CORRECTIONS ,    // Caluclate corrections, require dipole amplitude with rotational symmetry
    PRINT_NUCLEUS,
    F2,             // Calculate F2
    SATURATION_SCALE    // Print saturation scale on a g gripd
};

enum WAVEF
{
    GAUSLC,
    BOOSTEDGAUSSIAN,
    DVCS,
    NRQCD
};

std::string InfoStr(MODE mode);

vector<double> NRQCD_parameters_from_file(int id);

WAVEF wavef_model = BOOSTEDGAUSSIAN;

int MCpoints(double t);  // Automatic mc points


int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler); // Do not let code to crash if an error occurs, we should
    // have error handling everywhere!
    double maxr=99;    
    double Qsqr=0;
    double xbj=0; // x for F2
    double t=0.1;
    int A=1;
//    DGLAPDist *gd=0;  // Initialized and used if we have nucleus consisting of ipsatnucleons
    int he3_id=-1;   // Used to set He3 configuration
    double mint=0;
    double maxt=1.5;
    double tstep=0.1;
    std::vector<double> tlist;
    double maxb = 10. / 0.19733;    // GeV^-1
    int nbperp = 25;
    //double xpom=0.000959089;
    double w = 100;
    double xp = -1;
    bool skewedness = false;
    double qsfluct_sigma=0;
    Fluctuation_shape fluctshape = FLUCTUATE_QUARKS;
    bool auto_mcintpoints = false;
    std::string wavef_file = "";
    int rng_offset=0;
    double t_in_xpom=1.0;  // This multiplies t in the expression for xpom, if 0, then xpom is independent of t

    bool ipglasma=false;
    bool periodic_boundary_conditions=false;

    // nrqcd parameters
    double NRQCD_A=0.213;
    double NRQCD_B=-0.0157;
    int NRQCD_param_id = -1; // if >0, use specific parameters from datafile


    cout << "# SubNucleon Diffraction by H. Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2024" << endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl; 
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;

    if (string(argv[1])=="-help")
    {
        cout << "-Q2, -W, -xp: set kinematics" << endl;
        cout << "-dipole A [ipglasma,ipglasma_binary,ipsatproton,smoothnuke] [ipglasmafile ipglasmastep (fm), ipsat_proton_width ipsat_proton_quark_width] [fluxtube tube_normalization] [com]    com: move origin to Center of Mass (with constituent quark ipsat)" << endl;
        cout << "-corrections: calculate correction R_g^2(1+\beta^2) as a function of t. Requires rot. sym. dipole amplitude." << endl;
        cout << "-mcintpoints points/auto" << endl;
        cout << "-skewedness: enable skewedness in dipole amplitude" << endl;
        cout << "-qsfluct sigma: set width of Q_s fluctuations (0: disable); only for ipsatproton!" << endl;
        cout << "-qsfluctshape [local,quarks]: set Q_s^2 to fluctuate at each point / for each quark" << endl;
        cout << "-satscale: print saturation scale" << endl;
        cout << "-F2 Qsqr x: calculate structure function" << endl;
        cout << "-wavef_file filename" << endl;
        cout << "-wavef gauslc/boostedgaussian/DVCS/NRQCD" << endl;
        cout << "-He3 [config_id], REQUIRES A=3!"<< endl;
        cout << "-mint, -maxt, -tstep" << endl;
        cout << "-maxb, -nbperp" << endl;
        cout << "-nrqcd_parameters A B" << endl;
        cout << "-nrqcd_parameters_from_file" << endl;
        cout << "-periodic_boundary_conditions: use periodic boundary conditions" << endl;
        cout << "-no_t_in_xpom: do not include t dependence in xpom" << endl;
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
        else if (string(argv[i])=="-xp")
            xp=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-wavef")
        {
            if (string(argv[i+1])=="gauslc")
                wavef_model = GAUSLC;
            else if (string(argv[i+1])=="boostedgaussian")
                wavef_model = BOOSTEDGAUSSIAN;
            else if (string(argv[i+1])=="NRQCD")
                wavef_model = NRQCD;
            else if (string(argv[i+1])=="DVCS")
                wavef_model = DVCS;
            else
            {
                cerr << "Unknown wave function " << argv[i+1] << endl;
                exit(1);
                
            }
        }
        else if (string(argv[i])=="-dipole")
        {
            A = StrToInt(argv[i+1]);
            
            if (A==1)
            {
                if (string(argv[i+2])=="ipsatproton" or string(argv[i+2])=="ipsatprotonparam" or string(argv[i+2])=="lcpt")
                {
					if (string(argv[i+2])=="ipsatproton")
	                    amp = new Ipsat_Proton(MZSAT);
                    else if (string(argv[i+2])=="ipsatprotonparam")
                    {
                        IPsat_fit_parameteters ipsatparam;
                        ipsatparam.saturation=true;
                        // In this case the parameterers followinb Bp and BG are 
                        //    C, mu0, lambdag, Ag, mc
                        ipsatparam.C=StrToReal(argv[i+6]);
                        ipsatparam.mu0 = std::sqrt(1.1); // StrToReal(argv[i+6]);
                        ipsatparam.lambdag=StrToReal(argv[i+7]);
                        ipsatparam.Ag=StrToReal(argv[i+8]);
                        ipsatparam.mc=StrToReal(argv[i+5]);

                        amp = new Ipsat_Proton(MZSAT, ipsatparam);
                    }
					else
						amp = new Ipsat_Proton(LCPT);
                    ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+3]));
                    ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+4]));

                    ((Ipsat_Proton *)amp)->SetShape(GAUSSIAN);
                    if (argc > i + 5)
                    {
                        if (string(argv[i + 5]) == "fluxtube")
                        {
                            ((Ipsat_Proton *)amp)->SetStructure(CENTER_TUBES);
                            ((Ipsat_Proton *)amp)->SetFluxTubeNormalization(StrToReal(argv[i + 5]));
                        }
                        else if (string(argv[i + 5]) == "com")
                            ((Ipsat_Proton *)amp)->SetQuarkCenterOfMassToOrigin(true);
                    }
                }
                else if (string(argv[i+2])=="ipglasma")
                {
                    ipglasma=true;
                    amp = new IPGlasma(argv[i+3], StrToReal(argv[i+4]), TEXT);
                }
                else if (string(argv[i+2])=="ipglasma_binary")
                {
                    ipglasma=true;
                    amp = new IPGlasma(argv[i+3], StrToReal(argv[i+4]), BINARY);
                }
                else
                {
                    cerr << "Unknown dipole " << argv[i+1] << endl;
                    return -1;
                }
            }
            else    // A>1
            {
                if (string(argv[i+2])=="smoothnuke")
                {
                    amp = new Smooth_ws_nuke(A);
                }
                else
                {
                    // Construct nucleus
                    std::vector<DipoleAmplitude* > nucleons;
//                  NOTE: we need to initialize these separately (which is slow) if we want independent e-b-e fluct
//                 DipoleAmplitude* nucleon;
//                        nucleon = new Ipsat_Proton(MZSAT);
//                       nucleon->SetSkewedness(skewedness);

                    for (int j=0; j<A; j++)
                    {
                      DipoleAmplitude *nucleon; 
                        if (string(argv[i+2])=="ipsatproton" or string(argv[i+2])=="ipsatprotonparam")
                        {   
                            if (string(argv[i+2])=="ipsatproton")
                                nucleon = new Ipsat_Proton(MZSAT);
                            else if (string(argv[i+2])=="ipsatprotonparam")
                            {   
                                IPsat_fit_parameteters ipsatparam;
                                ipsatparam.saturation=true;
                                // In this case the parameterers followinb Bp and BG are 
                                //    C, mu0, lambdag, Ag, mc
                                ipsatparam.C=StrToReal(argv[i+6]);
                                ipsatparam.mu0 = std::sqrt(1.1); // StrToReal(argv[i+6]);
                                ipsatparam.lambdag=StrToReal(argv[i+7]);
                                ipsatparam.Ag=StrToReal(argv[i+8]);
                                ipsatparam.mc=StrToReal(argv[i+5]);

                                nucleon = new Ipsat_Proton(MZSAT, ipsatparam);
                            } 
                            ((Ipsat_Proton*)nucleon)->SetProtonWidth(StrToReal(argv[i+3]));
                            ((Ipsat_Proton*)nucleon)->SetQuarkWidth(StrToReal(argv[i+4]));

                            ((Ipsat_Proton *)nucleon)->SetShape(GAUSSIAN);
                            if (argc > i + 5)
                            {
                                if (string(argv[i + 5]) == "fluxtube")
                                {
                                    ((Ipsat_Proton *)nucleon)->SetStructure(CENTER_TUBES);
                                    ((Ipsat_Proton *)nucleon)->SetFluxTubeNormalization(StrToReal(argv[i + 5]));
                                }
                                else if (string(argv[i + 5]) == "com")
                                {
                                    ((Ipsat_Proton *)nucleon)->SetQuarkCenterOfMassToOrigin(true);
                                }

                                else if (string(argv[i + 5]).substr(0, 1) != "-" and string(argv[i + 2]) != "ipsatprotonparam")
                                {
                                    cerr << "Unknown ipsatproton option " << argv[i + 5] << endl;
                                    exit(1);
                                }
                            }
                        }
                        else if (string(argv[i+2])=="ipglasma_binary")
                            nucleon = new IPGlasma(argv[i+3], StrToReal(argv[i+4]), BINARY);
                        else if (string(argv[i+2])=="ipglasma")
                            nucleon = new IPGlasma(argv[i+3], StrToReal(argv[i+4]), TEXT);
                        nucleons.push_back(nucleon);
                    }
                    amp = new Nucleons(nucleons);
                    ((Nucleons*) amp)->SetDeuteronStructure(NUCLEONS);
                } // End construct nucleus
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
        else if (string(argv[i])=="-totalcrosssections")
            mode = TOTALCROSSSECTION;
        else if (string(argv[i])=="-qsfluct")
            qsfluct_sigma = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-qsfluctshape")
        {
            if (string(argv[i+1])=="quarks")
                fluctshape = FLUCTUATE_QUARKS;
            else
            {
                cerr << "Unknown fluctuation type " << argv[i+1] << endl;
                exit(1);
            }
        }
        else if (string(argv[i])=="-maxr")
            maxr = StrToReal(argv[i+1]);
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
        else if (string(argv[i])=="-He3")
            he3_id = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-rng_offset")
            rng_offset = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-mint")
            mint=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxt")
            maxt=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-tstep")
            tstep=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-tlist")
            tlist = StrToList(argv[i+1]);
        else if (string(argv[i])=="-maxb")
            maxb=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-nbperp")
            nbperp=StrToInt(argv[i+1]);
        else if (string(argv[i])=="-no_t_in_xpom")
            t_in_xpom = 0.0;
        else if (string(argv[i])=="-nrqcd_parameters")
        {
            NRQCD_A=StrToReal(argv[i+1]);
            NRQCD_B=StrToReal(argv[i+2]);
        }
        else if (string(argv[i])=="-nrqcd_parameters_from_file")
        {
            std::vector<double> params=NRQCD_parameters_from_file(StrToInt(argv[i+1]));
            NRQCD_A = params[0];
            NRQCD_B = params[1];
        }
        else if (string(argv[i])=="-periodic_boundary_conditions")
            periodic_boundary_conditions=true;
     else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
    }

    if (tlist.size() == 0) {
        for (double t = mint; t < maxt; t += tstep)
            tlist.push_back(t);
    }

    // Initialize global random number generator
    const gsl_rng_type * rngtype;
    gsl_rng_env_setup();
    long int seed =  gsl_rng_default_seed +rng_offset;
    cout << "# initializing rng seed to: " << seed << endl;
    rngtype = gsl_rng_default;
    global_rng = gsl_rng_alloc(rngtype);
    gsl_rng_set(global_rng, seed);

    if (ipglasma and periodic_boundary_conditions)
        ((IPGlasma*)amp)->SetPeriodicBoundaryConditions(periodic_boundary_conditions);
    if (!ipglasma and periodic_boundary_conditions)
    {   
        cerr << "Only IPGlasma dipoles support periodic boundary conditions! " << endl;
        exit(1);
    }


    WaveFunction *wavef;
    if (wavef_model == GAUSLC)
    {
        if (wavef_file == "") wavef_file = "gaus-lc.dat";
        wavef = new GausLC(wavef_file);
        cout << "# " << *(GausLC*)wavef << endl;
    }
    else if (wavef_model == BOOSTEDGAUSSIAN)
    {
        if (wavef_file == "") wavef_file = "gauss-boosted.dat";
        wavef = new BoostedGauss(wavef_file);
        cout << "# " << *(BoostedGauss*)wavef << endl;
    }
    else if (wavef_model == NRQCD)
    {
        wavef = new NRQCD_WF(NRQCD_A, NRQCD_B);
        FACTORIZE_ZINT=true;
        cout << "# " << *(NRQCD_WF*)wavef << endl;
    }
    else if (wavef_model == DVCS)
    {
        wavef = new DVCSPhoton;
        cout << "# " << *(DVCSPhoton*)wavef << endl;
    }

    amp->SetSkewedness(skewedness);
    if (qsfluct_sigma > 0)
    {
        if (A==1)
        {
        ((Ipsat_Proton*)amp)->SetQsFluctuation(qsfluct_sigma);
        ((Ipsat_Proton*)amp)->SetFluctuationShape(fluctshape);
        }
        else
        {
            std::vector<DipoleAmplitude*> nucleons = ((Nucleons*)amp)->GetNucleons();
            for (unsigned int i=0; i<nucleons.size(); i++)
            {
                ((Ipsat_Proton*)nucleons[i])->SetQsFluctuation(qsfluct_sigma);
                ((Ipsat_Proton*)nucleons[i])->SetFluctuationShape(fluctshape);
            }
        }
    }
    if (A==3)
    {
        ((Nucleons*)amp)->SetHeId(he3_id);
    }

    amp->InitializeTarget();

    Diffraction diff(*amp, *wavef);
    diff.SetMaxR(maxr*5.068);

    cout << "# " << InfoStr(mode) << endl;
    //cout << "# " << *wavef << endl;

    double mp = 0.938;
    double mjpsi = wavef->MesonMass();

    if (mode == PRINT_NUCLEUS)
    {
        if (ipglasma)
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
                
                WilsonLine wl =((IPGlasma*)amp)->GetWilsonLine(x,y);
                double tr = wl.Trace().real();
             
                cout << y/5.068 << " " << x/5.068 << " " << ((IPGlasma*)amp)->Amplitude(0.01, origin, p) << " " << ((IPGlasma*)amp)->AmplitudeImaginaryPart(0.01, origin, p) << " " << ((IPGlasma*)amp)->Amplitude(0.01, p, p) << " " << 1.0 - tr/3.0 <<endl;
            }
         cout << endl;
        }
        }
        else
        { 
        
        double origin[2]={0,0};
        double max = 25;
        double min = -25;
        double step = 0.1;
        cout << "# x y N(0,(x,y)) T(b) " << endl;
        for (double y=min+step/2; y < max-step/2; y+=step)
        {
            for (double x=min+step/2; x < max-step/2; x+=step)
            {
                double p[2] = {x,y};
                
                double density = 0;
                if (A == 1)
                    density =((Ipsat_Proton*)amp)->Density(Vec(x,y));
                else
                {
                    std::vector<DipoleAmplitude*> nucleons = ((Nucleons*)amp)->GetNucleons();
                    std::vector<Vec> positions =((Nucleons*)amp)->GetNucleonPositions();
                    for (unsigned int i=0; i<nucleons.size(); i++)
                    {
                        density += nucleons[i]->Density(Vec(x,y)-positions[i]);
                    }
                }
                cout << y/5.068 << " " << x/5.068 << " " << amp->Amplitude(0.001, origin, p) << " " << density << endl;
            }
            cout << endl;
        }
        }
        
         
        
         
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

        if (xp < 0)
            cout << "# Amplitude as a function of t, Q^2=" << Qsqr << ", W=" << w << endl;
        else
            cout << "# Amplitude as a function of t, Q^2=" << Qsqr << ", xp=" << xp << endl;
        cout << "# t  amplitude [GeV^-2] columns: transverse real, transverse imag, longitudinal real, longitudinal imag" << endl;

        for (auto t: tlist)
        {
            double xpom;
            if (xp < 0)
                xpom = (mjpsi*mjpsi+Qsqr+t_in_xpom*t)/(w*w+Qsqr-mp*mp);
            else
                xpom = xp;
            if (xpom > 0.04)
            {
                cerr << "xpom = " << xpom << ", can't do this!" << endl;
            }
            
            if(auto_mcintpoints)
                MCINTPOINTS = MCpoints(t);
            
            cout.precision(5);
            std::complex<double> trans = diff.ScatteringAmplitude(xpom, Qsqr, t, T);
            std::complex<double> lng(0.0,0.0);
            if (Qsqr > 0)
                lng = diff.ScatteringAmplitude(xpom, Qsqr, t, L);

            cout << t << " ";
            cout.precision(10);
            cout << trans.real()  << " " << trans.imag() << " " << lng.real() << " " << lng.imag() << endl;
        }
    }

    else if (mode == TOTALCROSSSECTION)
    {
        if (xp < 0)
            cout << "# Amplitude as a function of t, Q^2=" << Qsqr << ", W=" << w << endl;
        else
            cout << "# Amplitude as a function of t, Q^2=" << Qsqr << ", xp=" << xp << endl;
            
            double xpom;
            if (xp < 0)
            xpom = (mjpsi*mjpsi+Qsqr+t_in_xpom*t)/(w*w+Qsqr-mp*mp);
            else
            xpom = xp;
            if (xpom > 0.04)
            {
                cerr << "xpom = " << xpom << ", can't do this!" << endl;
            }
            
        auto data = diff.ComputeTotalCrossSection(xpom, Qsqr, nbperp, maxb);
        cout << "# Total cross section (transverse): " << data.sigma_T << " nb" << endl;
        cout << "# Total cross section (longitudinal): " << data.sigma_L << " nb" << endl;
        cout << "# b (GeV^-1)  F  columns: transverse real, transverse imag, longitudinal real, longitudinal imag, transverse |F|^2, longitudinal |F|^2" << endl;
        for (int ib=0; ib<nbperp; ++ib) {
            const double bval = data.b[ib];
            const std::complex<double>& FT = data.F_T[ib];
            std::complex<double> FL(0.,0.);
            double FT_sqr = data.F_T_sqr[ib];
            double FL_sqr = 0.0;
            if (Qsqr > 0) {
                FL = data.F_L[ib];
                FL_sqr = data.F_L_sqr[ib];
            }
            cout.precision(5);
            cout << std::fixed << bval << " ";
            cout.precision(10);
            cout << std::scientific
                 << FT.real() << " " << FT.imag() << " "
                 << FL.real() << " " << FL.imag() << " "
                 << FT_sqr << " " << FL_sqr << endl;
        }
    }

    else if (mode == CORRECTIONS)
    {
        cout << "# Real part correction" << endl;
        cout << "# t  transverse  longitudinal" << endl;
        double tstep=0.02;
        //for (t=mint; t<=maxt; t+=tstep)
        for (auto t: tlist)
        {
            double xpom = (mjpsi*mjpsi+Qsqr+t_in_xpom*t)/(w*w+Qsqr-mp*mp);
            if (xpom > 0.04)
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
        FACTORIZE_ZINT=true;
        cout << "#F2(Qsqr=" << Qsqr << ", xbj=" << xbj << "): light charm tot F_L(light) F_L(charm) F_L(tot)" << endl;
        double orig_x = xbj;
        WaveFunction * photon = new VirtualPhoton();;
        ((VirtualPhoton*)photon)->SetQuark(LIGHT, 0.03);
        cout << "# Quarks: " << ((VirtualPhoton*)photon)->GetParamString() << endl;
        
        amp->SetSkewedness(false);
        Diffraction f2(*amp, *photon);
       	f2.SetMaxR(maxr*5.068);
        cout << "#Maxr = " << f2.MaxR() << endl;
        // Use the fact that photon-proton cross section is just diffractive amplitude at t=0
        // Note* 4pi, as convention in BoostedGaussian and VirtualPhoton classes are different!!!
        double xs_t = f2.ScatteringAmplitude(xbj, Qsqr, 0, T).real();
        double xs_l = f2.ScatteringAmplitude(xbj, Qsqr, 0, L).real();
        double structurefun = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
        double fl_light =Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l;
        
        double mc=1.4;
        // heavy quark contribution
        ((VirtualPhoton*)photon)->SetQuark(C, mc);
        double xbj_c = xbj * (1.0 + 4.0*mc*mc / Qsqr);
        double xs_t_c = 0;
        double xs_l_c = 0;
        double fl_c = 0;
        double structurefun_c = 0;
        if (xbj_c < 0.01 or true)
        {
            cout << "# Quarks: " << ((VirtualPhoton*)photon)->GetParamString() << endl;
            xs_t_c = f2.ScatteringAmplitude(xbj_c, Qsqr, 0, T).real();
            xs_l_c = f2.ScatteringAmplitude(xbj_c, Qsqr, 0, L).real();
            structurefun_c = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_c+xs_t_c);
            fl_c =Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_c);
        }
        
        // b quark contribution
        ((VirtualPhoton*)photon)->SetQuark(B, 4.75);
        double xbj_b = xbj * (1.0 + 4.0*4.75*4.75 / Qsqr);
        double xs_t_b = 0;
        double xs_l_b = 0;
        double structurefun_b = 0;
        double fl_b = 0;
        if (xbj_b < 0.01 and false)
        {
            cout << "# Quarks: " << ((VirtualPhoton*)photon)->GetParamString() << endl;
            xs_t_b = f2.ScatteringAmplitude(xbj_b, Qsqr, 0, T).real();
            xs_l_b = f2.ScatteringAmplitude(xbj_b, Qsqr, 0, L).real();
            structurefun_b = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_b+xs_t_b);
            fl_b =Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_c);
        }
	
        cout << orig_x << " " << Qsqr << " " << structurefun << " " << structurefun_c << " " << structurefun + structurefun_c + structurefun_b << " " << fl_light << " " << fl_c << " " << fl_c + fl_light + fl_b<< endl;
        
        delete photon;
    }
    
    
    gsl_rng_free(global_rng);
    delete amp;
    delete wavef;
//    if (gd != 0)
//        delete gd;
}


string InfoStr(MODE mode)
{
    stringstream info;

    info << "Parameters: MCINTPOINTS: " << MCINTPOINTS << " ZINT_INTERVALS " << ZINT_INTERVALS << " MCINTACCURACY " << MCINTACCURACY << " ZINT_RELACCURACY " << ZINT_RELACCURACY;
    info << ". Integration method Suave ";

    info << endl << amp->InfoStr();

    if (FACTORIZE_ZINT)
        info <<"# z integral factorized";
    else info << "# z integral not factorized";

    return info.str();
}

int MCpoints(double t)
{
    if (t<0.1)
        return 5e5;
    else if (t<0.6)
        return 1e6;
    else
        return 1e7;
}


std::vector<double> NRQCD_parameters_from_file(int id)
{
    std::ifstream file("nrqcd_jpsi_parameters.dat");
    int i=0;
    std::vector<double> P;
    while(!file.eof())
    {
        string line;
        getline(file,line);
        if (line.substr(0,1)=="#")
            continue;
        if (i==id)
        {
            stringstream ss(line);
            double tmp; ss >>tmp;
            P.push_back(tmp);
            ss >> tmp;
            P.push_back(tmp);
            return P;
        }
        i=i+1;
    }
    std::cerr << "Unknown id " << id << " for NRQCD_parameters_from_file" << endl;
    exit(1);
    return P;
}

