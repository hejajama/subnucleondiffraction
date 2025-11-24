/* HERA sigmar fitter
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 *
 * Computes sigmar at given HERA kinematics
 *
 * 1st argument: JIMWLK solution dir
 * 2nd argument: step size
 * 3rd argument: max steps
 * 4th argument: HERA datafile
 *
 * HERA datafile with structure
 * Q^2 x inelasticity sigma_r error
 * Lines starting with # are ignored
 * inelasticity is the DIS y variable y=Q^2/(sx)
 *
 * Prints the computed reduced cross section values in the standard output, and the format is
 * Q^2 x y experimental_sigmar exp_error theory_sigmar
 */


#include "../src/dipole.hpp"
#include "../src/diffraction.hpp"


#include "../src/vector.hpp"
#include "../src/subnucleon_config.hpp"
#include "../src/ipglasma.hpp"
#include "../src/dis.hpp"
#include <amplitudelib/virtual_photon.hpp>
#include "../src/gitsha1.h"

#include <gsl/gsl_rng.h>

#include <tools/tools.hpp>
#include <gsl/gsl_errno.h>

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sstream>
#include <gsl/gsl_errno.h>
using namespace std;
using namespace Amplitude;

//double x0 = 0.01;
double x0 = 0.01; //0.00175;
double ds = 0.0004;
gsl_rng* global_rng;
WilsonLineDataFileType DATATYPE = BINARY;

double mean(vector<double> &v)
{
    double sum=0;
    for (unsigned int i=0; i<v.size(); i++) sum += v[i];
    return sum/v.size();
}

int main(int argc, char* argv[])
{
	double minq2=1;
    MCINTPOINTS=5e4; 
//    MCINTPOINTS=1e5;
    FACTORIZE_ZINT=true;
    gsl_set_error_handler(&ErrHandler);
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
	double ipglasma_step = 5.12/700.0; //0.01;
    cout << "# F2 fitter, intpoints " << MCINTPOINTS <<  endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl;
    if (argc < 8)
    {
        cout << "Syntax: " << argv[0] << " jimwlkdir step maxstep heradata alphas config ds x0 schwinger"  << endl;
        return 0;
    }
    
    string jimwlkdir = argv[1];
    int step = StrToInt(argv[2]);
    int maxstep = StrToInt(argv[3]);
    double alphas = StrToReal(argv[5]);
    string herafile = argv[4];
    double quarkmass = 1.4;
    int averages = StrToInt(argv[6]);
    ds = StrToReal(argv[7]); 
    x0  = StrToReal(argv[8]);
    double schwinger = StrToReal(argv[9]);
   
    cout << "# Command: " ; for (unsigned int i=0; i<argc; i++) cout << argv[i] << " ";
    cout << endl;
 
    bool fixed_coupling = false;
    double minx=0;
    if (alphas > 0)
    {
        // Fixed coupling
        fixed_coupling = true;
        double maxy =  M_PI*M_PI / alphas * maxstep * ds;
        minx = x0*exp(-maxy);
        cout << "# x range in fit: " << x0 << " - " << minx << ", # of configs" << averages << endl;
    }
    else
    {
        double maxy =  M_PI*M_PI * maxstep * ds;
        minx = x0*exp(-maxy);
        cout << "# x range in fit: " << x0 << " - " << minx << ", # of configs" << averages << endl;

    }
    
    bool scale_x = true;   // scale bjorken x to take into account the quark mass
    bool include_light = false;
    
    int points=0;
    double chisqr=0;
   
    cout << "# "; if (include_light) cout << "Light quarks included. "; cout << " Charm mass " << quarkmass << endl; 
    
    if (scale_x)
        cout <<"# Shifting Bjorken-x to x*(1 + 4mf^2/Q^2)" << endl;
    
    // Rear HERA data
    vector<double> xvals, yvals, qsqrvals, expvals, experrors;
    //NB: yvals is now vector of the inelasticities, y=Q^2/(sx)
    ifstream file(herafile.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << herafile.c_str() << endl;
        return -1;
    }
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#")
        continue;
        string x,qsqr,y,sigmar,err;
        stringstream l(line);
        l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
        //if (std::abs(StrToReal(f2)) <1e-10) continue;
        if (StrToReal(x)>x0 or StrToReal(x)<minx) continue;
        qsqrvals.push_back(StrToReal(qsqr)); xvals.push_back(StrToReal(x));
        yvals.push_back(StrToReal(y));
        expvals.push_back(StrToReal(sigmar)); experrors.push_back(StrToReal(err));
    }
    file.close();
   	cout << "# Loaded " << expvals.size() << " datapoints from file " << herafile << endl; 
    
    
    
    cout << "# Q^2 [GeV^2]  x  y  HERA-\\sigma_r  HERA-err theory-\\sigma_r " << endl;
    // Compute reduced cross section
    for (unsigned int i=0; i<xvals.size(); i++)
    {
	if (qsqrvals[i] < minq2)
		continue;
        double x = xvals[i];
        double sqrts = std::sqrt( qsqrvals[i]/(x * yvals[i]) );
      	double xcharm = x; 
        if (scale_x)
            xcharm = x   * (1.0 + 4.0*SQR(quarkmass)/qsqrvals[i]);
       
		if (x > x0 or xcharm > x0)
			continue;
        // Calculate reduced xs below and above, and interpolate
        // Fixed coupling: step = alphas*ln(x0/x) / (pi^2 * ds)
        int steps_upper = 0;
        int steps_lower = 0;
        double evolsteps;

		// Charm
		int steps_upper_c = 0;
		int steps_lower_c = 0;
		double evolsteps_c;
         
        if (fixed_coupling)
        {
            evolsteps = alphas * log(x0/x)/ (M_PI*M_PI*ds);
	    evolsteps_c = alphas * log(x0/xcharm)/ (M_PI*M_PI*ds);
		
	}
	else 
	{
		evolsteps =  log(x0/x)/ (M_PI*M_PI*ds);
                evolsteps_c =  log(x0/xcharm)/ (M_PI*M_PI*ds);
	}

	int evolstepsint = (int)(evolsteps );
	int evolstepsint_c = (int)(evolsteps_c );
			
            steps_lower = evolstepsint - (evolstepsint%step) ;
            steps_upper = evolstepsint + (step - evolstepsint%step) ;

			steps_lower_c = evolstepsint_c - (evolstepsint_c%step);
			steps_upper_c = evolstepsint_c + (step - evolstepsint_c%step); 
            //cout << " Evolsteps " << evolsteps << " lower " << steps_lower << " upper " << steps_upper << endl;
	cout << "# Evolsteps_c " << evolsteps_c << " lower " << steps_lower_c << " upper " << steps_upper_c << endl;
        
        VirtualPhoton charm;
        charm.SetQuark(Amplitude::C, quarkmass);
		VirtualPhoton light;
        light.SetQuark(Amplitude::LIGHT, 0.14);
        vector<double> xs_t_upper;
        vector<double> xs_l_upper;
        vector<double> xs_t_lower;
        vector<double> xs_l_lower;

		vector<double> xs_t_upper_c;
        vector<double> xs_l_upper_c;
        vector<double> xs_t_lower_c;
        vector<double> xs_l_lower_c;


        //for(int conf = 0; conf < averages; conf++)
        int conf = averages; // Do only single configuration, average later separately!
        {
            stringstream fname_upper;
            stringstream fname_lower;
            stringstream fname_upper_c;
            stringstream fname_lower_c;
            if (DATATYPE == TEXT)
            {
                fname_upper << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_upper;
                fname_lower << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_lower;
			
                fname_upper_c << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_upper_c;
                fname_lower_c << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_lower_c;
             } else if (DATATYPE == BINARY)
             {
                fname_upper << jimwlkdir << "/V-" << conf << "_steps_" << steps_upper;
                fname_lower << jimwlkdir << "/V-" << conf << "_steps_" << steps_lower;
			
                fname_upper_c << jimwlkdir << "/V-" << conf << "_steps_" << steps_upper_c;
                fname_lower_c << jimwlkdir << "/V-" << conf << "_steps_" << steps_lower_c;

             }	



            //cout << fname_upper_c.str() << " " << fname_lower_c.str() << endl;
            
	  
		     
 			if (include_light)
			{
	            IPGlasma dipole_upper(fname_upper.str(), ipglasma_step, DATATYPE);
    	        IPGlasma dipole_lower(fname_lower.str(), ipglasma_step, DATATYPE);


				if (schwinger > 0)
				{
					dipole_upper.SetSchwinger(true, schwinger);
					dipole_lower.SetSchwinger(true, schwinger);
				}


            	Diffraction f2upper(dipole_upper, light);
            	Diffraction f2lower(dipole_lower, light);
 
               	double t_up = f2upper.ScatteringAmplitude(x, qsqrvals[i], 0, T).real();
            	double l_up = f2upper.ScatteringAmplitude(x, qsqrvals[i], 0, L).real();
            	xs_t_upper.push_back(t_up);
            	xs_l_upper.push_back(l_up);
            
                double t_low = f2lower.ScatteringAmplitude(x, qsqrvals[i], 0, T).real();
                double l_low = f2lower.ScatteringAmplitude(x, qsqrvals[i], 0, L).real();
            	xs_t_lower.push_back(t_low);
            	xs_l_lower.push_back(l_low);
            }
			else
			{
				xs_t_upper.push_back(0); xs_l_upper.push_back(0);
				xs_t_lower.push_back(0); xs_l_lower.push_back(0);
			}
          	
		IPGlasma dipole_upper_c(fname_upper_c.str(), ipglasma_step, DATATYPE);
            IPGlasma dipole_lower_c(fname_lower_c.str(), ipglasma_step, DATATYPE);

            if (schwinger > 0)
            {
				dipole_upper_c.SetSchwinger(true, schwinger);
                dipole_lower_c.SetSchwinger(true, schwinger);

            }
         
			Diffraction f2upper_c(dipole_upper_c, charm);
            Diffraction f2lower_c(dipole_lower_c, charm);
            

            double t_up_c = f2upper_c.ScatteringAmplitude(xcharm, qsqrvals[i], 0, T).real();
            double l_up_c = f2upper_c.ScatteringAmplitude(xcharm, qsqrvals[i], 0, L).real();
            xs_t_upper_c.push_back(t_up_c);
            xs_l_upper_c.push_back(l_up_c);

            
            double t_low_c = f2lower_c.ScatteringAmplitude(xcharm, qsqrvals[i], 0, T).real();
            double l_low_c = f2lower_c.ScatteringAmplitude(xcharm, qsqrvals[i], 0, L).real();
            xs_t_lower_c.push_back(t_low_c);
            xs_l_lower_c.push_back(l_low_c);

 
            
        }
        
        double xs_t_u =mean(xs_t_upper);
        double xs_l_u = mean(xs_l_upper);
        double f2_upper = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_u + xs_l_u);
        double fl_upper = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_u );
        double sigmar_upper = f2_upper - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_upper;

        double xs_t_l =mean(xs_t_lower);
        double xs_l_l = mean(xs_l_lower);
        double f2_lower = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_l + xs_l_l);
        double fl_lower = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_l );
        double sigmar_lower = f2_lower - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_lower;
        
        // Interpolate in s ~ ln 1/x
        double sigmar = sigmar_lower + (evolsteps-(double)(steps_lower)) / ((double)(steps_upper - steps_lower)) * (sigmar_upper - sigmar_lower);
        

		double xs_t_u_c =mean(xs_t_upper_c);
        double xs_l_u_c = mean(xs_l_upper_c);
        double f2_upper_c = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_u_c + xs_l_u_c);
        double fl_upper_c = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_u_c );
        double sigmar_upper_c = f2_upper_c - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_upper_c;
       
	    
        double xs_t_l_c =mean(xs_t_lower_c);
        double xs_l_l_c = mean(xs_l_lower_c);
        double f2_lower_c = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_l_c + xs_l_l_c);
        double fl_lower_c = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_l_c );
        double sigmar_lower_c = f2_lower_c - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_lower_c;
       
	   	
	  cout << "#f2_c_upper " << f2_upper_c << " f2_lower_c " << f2_lower_c << endl;
	cout << "#f2_upper " << f2_upper << " f2_lower " << f2_lower << endl;

	    
        // Interpolate in s ~ ln 1/x
        double sigmar_c = sigmar_lower_c + (evolsteps_c-(double)(steps_lower_c)) / ((double)(steps_upper_c - steps_lower_c)) * (sigmar_upper_c - sigmar_lower_c);
		cout << "#sigmar_lower " << sigmar_lower_c << " upper " << sigmar_upper_c << endl;        

	double f2_light = f2_lower + (evolsteps-(double)(steps_lower)) / ((double)(steps_upper - steps_lower)) * (f2_upper - f2_lower);
	double f2_c = f2_lower_c + (evolsteps_c-(double)(steps_lower_c)) / ((double)(steps_upper_c - steps_lower_c)) * (f2_upper_c - f2_lower_c);

        cout << qsqrvals[i] << " " << xvals[i] << " " << yvals[i] << " " << expvals[i] << " " << experrors[i] << " " << sigmar + sigmar_c  << endl; //<<  " " << f2_light << " " << f2_c << endl << "#" << endl;

        points++;
        chisqr += pow((sigmar - expvals[i])/experrors[i],2.0);
    }
    
    cout << "#Chi^2/N " << chisqr/points << " points " << points << endl;
    
}
