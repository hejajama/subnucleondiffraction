/* Calculate avarage dipole amplitude averaged over b
*/

#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include "../src/subnucleon_config.hpp"
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <gsl/gsl_rng.h>

gsl_rng* global_rng;
using namespace std;

struct bhelper
{
    IPGlasma *glasma;
    double r;
    double theta_b;
    double b;
    double theta_r_b;
};

const int INTPOINTS = 6;
const double INTACCURACY = 0.1;
double B_BIN=0.1;
double rhelperf(double theta_r, void* p); // average over dipole orientation
double bhelperf(double b, void* p); // average over b
double bhelperf_theta(double theta_b, void* p);

int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    double r= StrToReal(argv[3]);
    double b = StrToReal(argv[2]);
    double step = StrToReal(argv[4]);
    cout << "# Filename: " << fname << " b " << b << " Gev^-1 r " << r << "  step[fm]" << step <<  endl;

    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
   

    // If file extension is txt, 
    WilsonLineDataFileType wlinetype;
    if(fname.substr( fname.length() - 4 ) == ".txt") 
        wlinetype=TEXT;
    else
        wlinetype=BINARY;
    IPGlasma glasma(fname, step, wlinetype);
    bhelper helper;
    helper.glasma = &glasma;
	helper.r=r;
    
    gsl_function f;
    f.function = &bhelperf;  // do a circle around the proton
	f.params = &helper;
    
   cout << "#theta_r_b   N" << endl; 
    for (double theta_r_b=0; theta_r_b <= 2.0*M_PI; theta_r_b += 2.0*M_PI/30.0)
    {
    
        helper.theta_r_b = theta_r_b;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
	//MAXB = (5.12/2.0*5.068-r/2.0);	// so that q/qbar is never outside the lattice
        double result,error;
        int status = gsl_integration_qag(&f, b, b+B_BIN, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
	//int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 1e-5, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        //if (status)
          //  cerr << "#bint failed (" << gsl_strerror (status) << ", result  " << result << " relerror " << error/result << " r " <<r  << endl;
        
        
        gsl_integration_workspace_free(w);
  		//result /= (2.0*M_PI);     
        result = result / (M_PI*( (b+B_BIN)*(b+B_BIN) - b*b));
	//	result = result / (2.0*M_PI*2.0*M_PI);
        
        cout << theta_r_b << " " << result << endl;
    }
    
    


}


double bhelperf(double b, void* p)
{
    bhelper* par = (bhelper*)p;
    par->b = b;
    gsl_function f;
    f.function = &bhelperf_theta; // do a circle around the proton
    f.params = par;
    
    
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    //if (status)
      //  cerr << "#btheta failed (" << gsl_strerror (status)  << "), result  " << result << " relerror " << error/result << " b " <<b << endl;
    
    gsl_integration_workspace_free(w);
    
    return b*result;
    
    
}


double bhelperf_theta(double theta_b, void* p)
{
    bhelper* par = (bhelper*)p;
    par->theta_b = theta_b;
    gsl_function f;
    f.function = &rhelperf;  // rotate dipole
    f.params = par;
   
    double theta_r = theta_b - par->theta_r_b;  // theta_b_r is theta_b - theta_r 
	return rhelperf(theta_r, par); 

    
    

    
}

double rhelperf(double theta_r, void* p)
{
    bhelper* par = (bhelper*)p;
    
    Vec b(par->b * cos(par->theta_b), par->b * sin(par->theta_b));
    Vec r(par->r * cos(theta_r), par->r * sin(theta_r));
    Vec r2 = r*0.5;
    //Vec q1 = b; // + r2;
    //Vec q2 = b + r; //- r2;
    Vec q1 = b + r2;
    Vec q2 = b - r2;
   
    return par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
}
