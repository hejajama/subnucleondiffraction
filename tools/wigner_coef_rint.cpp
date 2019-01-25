/* Calculate rint dphi int for Wigner coefficients https://arxiv.org/pdf/1609.05773.pdf eqs 19,20
 */ 
#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include <tools/tools.hpp>
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <tools/tools.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

gsl_rng* global_rng;
using namespace std;

struct bhelper
{
    IPGlasma *glasma;
    double r;
    double theta_b;
    double b;
    double theta_r_b;
    double k;
};

const int INTPOINTS = 6;
const double INTACCURACY = 0.0001;
double helperf_theta(double theta_b, void* p);
double helperf_r(double r, void* p);

double inthelperf_mc( double *vec, size_t dim, void* par);

int MCPOINTS = 1e6;
bool MONTECARLO = true;;

int wigner_coef = 0;     //1=xW1, 0 = xW0

gsl_integration_workspace *w_overall;
gsl_integration_workspace *w_relative;


int main(int argc, char* argv[])
{
    const double fmgev=5.068;
    // Arguments: ipglasma filename step  cos2phimom
    if (argc < 3)
    {
        cerr << "Arguments: filaneme wlinestepsize [fm] " << endl;
        return 1;
    }
    string fname = argv[1];
    double step = StrToReal(argv[2]);
    
    if (MONTECARLO)
        {
	MCPOINTS = (int)(StrToReal(argv[3]));
    	cout << "# MC, intpoints=" << MCPOINTS << endl;
}

    cout << "# Filename: " << fname <<  "  step[fm] " << step <<  endl;
    double step_gev = step *fmgev;
/*    double maxb = 1.11*fmgev;
    double bstep = 0.05 * fmgev;
    double maxk = 4.01;
    double kstep = 0.25;
  */

/*    double minb = 0;
    double maxb = 0.61*fmgev;
    double bstep = 0.02 * fmgev;
 */
    double minb = 0.325*fmgev;
    double maxb = 0.475*fmgev+0.01;  
    double bstep = 0.025*fmgev;

/* Large k
    double mink = 0.2609196; 
    double maxk = 15.01;
    double kstep = 1.2;
*/

/* Smallish k*/
double mink = 0.01; //0.217433*1.2; // 1.6155459071999998; //0.217433;
double maxk = 13.50;
double kstep = 1.2;
// Linear large-k -part
/*
double mink = 0.25;
double kstep = 0.25;
double kmax = 3.1;
double bstep = 0.02 * fmgev;
double minb = 0.34*fmgev;
double maxb = 0.461*fmgev;
*/
   // start
/*
   double mink = 0.01;
   double kstep_start = 1.6;
   double kstep2 = 1.2;
   double kmax = 0.3;
*/

// continue
/*
    double mink=0.3131028;
       double kstep_start=1.6;
       double kstep2 = 1.3;
       double kmax=14.9;
*/
	/*
	mink=1.25;
	kstep=0.25;
	kmax=3.1;
     */
  
/*	mink=0;
	kstep=0.25;
	kmax=3.1;
*/
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    WilsonLineDataFileType wlinetype;
    if(fname.substr( fname.length() - 4 ) == ".txt")
    wlinetype=TEXT;
    else
    wlinetype=BINARY;
    
    IPGlasma glasma(fname, step, wlinetype);
    bhelper helper;
    helper.glasma = &glasma;
    
    gsl_function f;
    f.function = &helperf_r;  // do a circle around the proton
    f.params = &helper;
    
    w_relative = gsl_integration_workspace_alloc(INTPOINTS);
    w_overall = gsl_integration_workspace_alloc(INTPOINTS);
    
    // MC
    double *lower = new double[3];
    double *upper = new double[3];
    
    
    lower[0]=lower[1]=lower[2]=0;
    upper[0] = 5*5.068 ; // Max r
    upper[1] = 2.0*M_PI; // theta_r_angle
    upper[2] = 2.0*M_PI; // overall rotation
    
    gsl_monte_function mc_f;
    mc_f.f = &inthelperf_mc;
    mc_f.dim = 3;
    
    mc_f.params = &helper;
    mc_f.params = &helper;
    
    
    
    
    
    
    cout << "#k  b   integral" << endl;
//    for (double k=mink; k <=kmax;  k+=kstep)
//    double kstep = kstep_start;
    for (double k=mink; k<=maxk; k*=kstep)
    {
//        if (k >= 0.1)
//            kstep = kstep2;
        for (double b=minb; b <= maxb; b+=bstep)
        {
            
            helper.k=k;
            helper.b=b;
            
            double result,error;
            
            if (MONTECARLO)
            {
                gsl_monte_miser_state *s = gsl_monte_miser_alloc(mc_f.dim);
                gsl_monte_miser_integrate(&mc_f, lower, upper, mc_f.dim, MCPOINTS, global_rng, s, &result, &error);
                gsl_monte_miser_free(s);
            }
            else
            {
                gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
                
                int status = gsl_integration_qag(&f, 0, 4*fmgev, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
                
                
                gsl_integration_workspace_free(w);
            }
            
            
            result /= (2.0*M_PI);
            cout << k << " " << b << " " << result << " " << error/(2.0*M_PI) << endl;
        }
    }
    
    
    return 0;
    
}

double helperf_r(double r, void* p)
{
    bhelper *par = (bhelper*)p;
    par->r = r;
    
    double result,error;
    gsl_function fun;
    fun.params=par;
    fun.function = helperf_theta;
    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w_overall, &result, &error);
    
    if (wigner_coef == 0)
    {
        result *= r * gsl_sf_bessel_J0(par->k*r);
    }
    else
    {
        result *= r * gsl_sf_bessel_Jn(2,par->k*r);
    }

    return result;
}

    
double helperf_rb(double theta_rb, void* p);
// overall angle
double helperf_theta(double theta_b, void* p)
{
    bhelper *par = (bhelper*)p;
    par->theta_b = theta_b;
    gsl_function f;
    f.params = par;
    f.function = &helperf_rb;
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w_relative, &result, &error);
    
   
    return result;
}
    
double helperf_rb(double theta_r_b, void* p)
{
        
    bhelper* par = (bhelper*)p;
    double theta_b = par->theta_b;
        
    Vec b(par->b * cos(theta_b), par->b * sin(theta_b));
    Vec r(par->r * cos(theta_b+theta_r_b), par->r * sin(theta_b + theta_r_b));
        
    Vec r2 = r*0.5;

    Vec q1 = b + r2;
    Vec q2 = b - r2;
    
    if (wigner_coef == 1)
    {
        
        return std::cos(2.0*theta_r_b) * par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
    }
    else
    {
        return par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
    }
        
        
}



double inthelperf_mc( double *vec, size_t dim, void* p)
{
    double r = vec[0];
    double theta_r_b = vec[1];
    double overall_rotation = vec[2];
    
    bhelper *par = (bhelper*) p;
    par->r = r;
    par->theta_b = overall_rotation;
    par->theta_r_b = theta_r_b;

      double result = 0;
    if (wigner_coef == 0)
    {
         result = r * gsl_sf_bessel_J0(par->k*r);
     }
     else
     {
         result = r * gsl_sf_bessel_Jn(2,par->k*r);
 }

    return result * helperf_rb(theta_r_b, par);
 
}

