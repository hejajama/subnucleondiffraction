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

#include <tools/tools.hpp>
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
#include <amplitudelib/virtual_photon.hpp>
#include "gitsha1.h"

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

enum WAVEF
{
    GAUSLC,
    BOOSTEDGAUSSIAN,
    NRQCD
};

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
    
    //double xpom=0.000959089;
    double w = 100;
    double xp = -1;
    bool skewedness = false;
    double qsfluct_sigma=0;
    Fluctuation_shape fluctshape = LOCAL_FLUCTUATIONS;
    bool auto_mcintpoints = false;
    std::string wavef_file = "";
    bool schwinger = false;
    double schwinger_rc = 0;
    int rng_offset=0;
    double t_in_xpom=0.0;  // This multiplies t in the expression for xpom, if 0, then xpom is independent of t
    double jimwlk_steps=-1; // converted to W
    
    bool ipglasma=false;
    bool periodic_boundary_conditions=false;
    
    // nrqcd parameters
    double NRQCD_A=0.213;
    double NRQCD_B=-0.0157;
    int NRQCD_param_id = -1; // if >0, use specific parameters from datafile
    
    
    cout << "# SubNucleon Diffraction by H. Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2021" << endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl; 
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    
    if (string(argv[1])=="-help")
    {
        cout << "-Q2, -W, -xp: set kinematics" << endl;
        cout << "-real, -imag: set real/imaginary part" << endl;
        cout << "-dipole A [ipglasma,ipglasma_binary,ipsatproton,smoothnuke] [ipglasmafile ipglasmastep (fm), ipsat_proton_width ipsat_proton_quark_width] [fluxtube tunbe_normalization] [albacete] [com]    com: move origi nto Center of Mass (with constituent quark ipsat)" << endl;
        cout << "-corrections: calculate correction R_g^2(1+\beta^2) as a function of t. Requires rot. sym. dipole amplitude." << endl;
        cout << "-mcintpoints points/auto" << endl;
        cout << "-skewedness: enable skewedness in dipole amplitude" << endl;
        cout << "-qsfluct sigma: set width of Q_s fluctuations (0: disable); only for ipsatproton!" << endl;
        cout << "-qsfluctshape [local,quarks]: set Q_s^2 to fluctuate at each point / for each quark" << endl;
        cout << "-satscale: print saturation scale" << endl;
        cout << "-F2 Qsqr x: calculate structure function" << endl;
        cout << "-wavef_file filename" << endl;
        cout << "-wavef gauslc/boostedgaussian/NRQCD" << endl;
        cout << "-He3 [config_id], REQUIRES A=3!"<< endl;
        cout << "-schwinger r_c" << endl;
        cout << "-mint, -maxt, -tstep" << endl;
        cout << "-nrqcd_parameters A B" << endl;
        cout << "-nrqcd_parameters_from_file" << endl;
        cout << "-periodic_boundary_conditions: use periodic boundary conditions" << endl;
        cout << "-mcint [miser,vegas]: slect MC integral algorithm" << endl;
        cout << "-t_in_xpom: do not include t dependence in xpom" << endl;
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
                if (string(argv[i+2])=="ipsatproton" or string(argv[i+2])=="lcpt")
                {
					if (string(argv[i+2])=="ipsatproton")
	                    amp = new Ipsat_Proton(MZSAT);
					else
						amp = new Ipsat_Proton(LCPT);
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
                            else if (string(argv[i+5])=="com")
                                 ((Ipsat_Proton*)amp)->SetQuarkCenterOfMassToOrigin(true);
                            else if (string(argv[i+5]).substr(0,1)!="-")
                            {
                                cerr << "Unknown ipsatproton option " << argv[i+4] << endl;
                                exit(1);
                            }
                        }
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
                        if (string(argv[i+2])=="ipsatproton")
                        {
                            //if (j==0)
                            //    gd = new DGLAPDist;
                            //Ipsat_Proton *nucleon = new Ipsat_Proton(gd);
                            nucleon = new Ipsat_Proton(MZSAT);
                            ((Ipsat_Proton*)nucleon)->SetProtonWidth(StrToReal(argv[i+3]));
                            ((Ipsat_Proton*)nucleon)->SetQuarkWidth(StrToReal(argv[i+4]));

   
                            if (string(argv[i+5]) == "ALBACETE")
                                ((Ipsat_Proton*)nucleon)->SetShape(ALBACETE);
                            else
                            {
                                ((Ipsat_Proton*)nucleon)->SetShape(GAUSSIAN);
                                if (argc > i+5)
                                {
                                    if (string(argv[i+5])=="fluxtube")
                                    {
                                        ((Ipsat_Proton*)nucleon)->SetStructure(CENTER_TUBES);
                                        ((Ipsat_Proton*)nucleon)->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                                    }
                                    else if (string(argv[i+5])=="com")
                                    {
                                         ((Ipsat_Proton*)nucleon)->SetQuarkCenterOfMassToOrigin(true); 
                                    }

                                    else if (string(argv[i+5]).substr(0,1)!="-")
                                    {
                                        cerr << "Unknown ipsatproton option " << argv[i+5] << endl;
                                        exit(1);
                                    }
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
        else if (string(argv[i])=="-schwinger")
        {
            schwinger = true; 
            schwinger_rc = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-rng_offset")
            rng_offset = StrToInt(argv[i+1]);
        else if (string(argv[i])=="-mint")
            mint=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxt")
            maxt=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-tstep")
            tstep=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-t_in_xpom")
            t_in_xpom = 1.0;
        else if (string(argv[i])=="-nrqcd_parameters")
        {
            NRQCD_A=StrToReal(argv[i+1]);
            NRQCD_B=StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-nrqcd_parameters_from_file")
        {
            std::vector<double> params=NRQCD_parameters_from_file(StrToInt(argv[i+1]));
            NRQCD_A = params[0];
            NRQCD_B = params[1];
        }
        else if (string(argv[i])=="-periodic_boundary_conditions")
            periodic_boundary_conditions=true;
        else if (string(argv[i])=="-mcint")
        {
            if (string(argv[i+1])=="miser")
                MCINT = MISER;
            else if (string(argv[i+1])=="vegas")
                MCINT = VEGAS;
            else
            {
                cerr << "Unknown MC algorithm " << argv[i+1] << endl;
                exit(1);
            }

        }
        else if (string(argv[i])=="-lhc")
            KINEMATICS = LHC;
        else if (string(argv[i])=="-rhic")
            KINEMATICS = RHIC;
        
        else if (string(argv[i])=="-ff")
        {
            if (string(argv[i+1])=="starlight")
                NUCLEAR_FF = STARLIGHT;
            else if (string(argv[i+1])=="pointcharge")
                NUCLEAR_FF = POINT_CHARGE;
            else
            {
                cerr << "Unknown form factor " << argv[i+1] << endl;
                exit(1);
            }
        }
        else if (string(argv[i])=="-jimwlk_steps")
        {
            jimwlk_steps = StrToReal(argv[i+1]);
            const double ds = 0.004;
            double tmpy = jimwlk_steps*ds*M_PI*M_PI;
            double tmp_xp = 0.01 * std::exp(-tmpy);
            w = std::sqrt(3.097*3.097/tmp_xp);

            cout << "#JIMWLK_step sets w=" << w << endl;
            if (w > 29 and w < 32) w = 31.52590399193987;
            if (w > 490 and w < 500) w = 493.1481109621738; 
            cout << "#Final W to be used: "<< w << endl;

        }
     else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
        
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
    

    if (schwinger) ((IPGlasma*)amp)->SetSchwinger(true, schwinger_rc);
    

    amp->InitializeTarget();
    

    Diffraction diff(*amp, *wavef);
    diff.SetMaxR(maxr*5.068);

    
    cout << "# " << InfoStr() << endl;
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
                
                WilsonLine &wl =((IPGlasma*)amp)->GetWilsonLine(x,y);
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
        
        cout << "# Amplitude, t=" << mint <<", Q^2=" << Qsqr << ", W=" << w << endl;
        cout << "# B M0 M1x M1y (re and im parts, note that M0 has different units!)" << endl;
        double t = mint;
        
        const int THPOINTS=6;
        const double MAXTH = M_PI;
        double bstep = 2.5; //5./3.;
        
        double MAXB;
        double MINB;
        if (KINEMATICS == LHC)
        {
            MAXB = 7000;
            MINB=2.0*6.62*5.068*0.9;


            //MINB = 7010;
            MAXB = 11000;
        }
       else
       {
           MAXB = 1000;
           MINB = 2.0*6.37*5.068*0.9;
       }
        
        const double IGNORE_THB_AFTER_B = 200;
        
        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        
        // LHC
        for (double B = MINB;  B < MAXB; B+=bstep)
        {
                if (B > 100) bstep = 10./1.5;
                if (B > 300) bstep = 50./1.5;
            //    if (KINEMATICS == RHIC and B > 500) bstep = 50;
                if (B > 1000) bstep = 100/1.5;

 
                
            
//#pragma omp parallel for
 //           for (int i=0; i<THPOINTS; i++)
 //           {
            for (double theta_B = 0; theta_B <= MAXTH; theta_B += MAXTH/THPOINTS) {

                //double theta_B = MAXTH/THPOINTS; * i;
                //double theta_B=1;
                //double theta_B = 2.0*M_PI * i;
                double xpom;
                if (xp < 0)
                    xpom = (mjpsi*mjpsi+Qsqr+t_in_xpom*t)/(w*w+Qsqr-mp*mp);
                else
                    xpom = xp;
                if (xpom > 0.04)
                {
                    cerr << "xpom = " << xpom << ", can't do this!" << endl;
                    //continue;
                }
                
                if(auto_mcintpoints)
                    MCINTPOINTS = MCpoints(t);
                
                cout.precision(10);
                std::vector<COMPONENT> complist;
                complist.push_back(M0);
                complist.push_back(M1x);
                complist.push_back(M1y);
                
                std::vector<double> results = diff.ScatteringAmplitude(xpom, Qsqr, t, B, theta_B, true, complist, T);
                std::vector<double> results_im = diff.ScatteringAmplitude(xpom, Qsqr, t, B, theta_B, false, complist, T);

                cout << B << " " << theta_B << " ";
                
                
                for (int i=0; i < results.size(); i++)
                {
                    cout << results[i] << " " << results_im[i] << " ";
                }
                cout << endl;
                
                // At large enough B we do not have th_B dependence anymore
                if (B > IGNORE_THB_AFTER_B)
                {
                    
                    for (theta_B = theta_B+MAXTH/THPOINTS; theta_B <= MAXTH; theta_B += MAXTH/THPOINTS)
                    {
                        cout << B << " " << theta_B << " ";
                        
                        
                        for (int i=0; i < results.size(); i++)
                        {
                            cout << results[i] << " " << results_im[i] << " ";
                        }
                        cout << endl;
                        
                    }
                        
                }
        }
            /*

                
                M0data[i] = results[0];
                M1xdata[i] = results[1];
                M1ydata[i] = results[2];
                
                M0data_im[i] = results_im[0];
                M1xdata_im[i] = results_im[1];
                M1ydata_im[i] = results_im[2];
             
                
            }
            
            M0data[THPOINTS] = M0data[0];
            M1xdata[THPOINTS] = M1xdata[0];
            M1ydata[THPOINTS] = M1ydata[0];
            
            M0data_im[THPOINTS] = M0data_im[0];
            M1xdata_im[THPOINTS] = M1xdata_im[0];
            M1ydata_im[THPOINTS] = M1ydata_im[0];
            
            
            
            std::cout << std::scientific;
            std::cout << std::setprecision(10);
            for (int i=0; i<THPOINTS+1; i++)
            {
                double theta_B = MAXTH/THPOINTS * i;
                cout << B << " ";
                cout << theta_B  << " ";
                cout << M0data[i]  << " " << M0data_im[i]
                    << " " << M1xdata[i] << " " << M1xdata_im[i]
                    << " "<< M1ydata[i] << " " << M1ydata_im[i] << endl;
                
            }
             */
           // exit(1);
                
        }
    }
    else if (mode == CORRECTIONS)
    {
        cout << "# Real part correction" << endl;
        cout << "# t  transverse  longitudinal" << endl;
        double tstep=0.02;
        for (t=mint; t<=maxt; t+=tstep)
        {
            double xpom = (mjpsi*mjpsi)/(w*w); //+Qsqr+t_in_xpom*t)/(w*w+Qsqr-mp*mp);
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
    
    
    
    gsl_rng_free(global_rng);


    
    delete amp;
    delete wavef;
    
//    if (gd != 0)
//        delete gd;
    
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
    
    info << endl << amp->InfoStr() << endl;;

    
    
    
    if (KINEMATICS == RHIC)
        info << "# RHIC kinematics";
    else if (KINEMATICS == LHC)
        info << "# LHC kinematics";
    info << endl;
    if (NUCLEAR_FF == STARLIGHT)
        info << "# STARLIGHT form factor";
    else if (NUCLEAR_FF == POINT_CHARGE)
        info << "# Point charge form factor";
    else if (NUCLEAR_FF == WOODSAXON)
        info << "# Woods Saxon form factor";
    else
        info << "# Unknown form factor ";
    
    info << endl << "# Off forward phase: " << OFF_FORWARD_PHASE << endl;
    
    return info.str();

}

int MCpoints(double t)
{
    if (t<0.1)
        return 5e6;
    else if (t<0.6)
        return 1e7;
    else
        return 1e8;
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
