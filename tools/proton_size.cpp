/*
 * Calculate proton size based on dipole amplitude
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

// Implement definition from 1407.8458 Eq. 17
// N(b=R_eff, r~r_0) = N0

#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include <tools/tools.hpp>
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <tools/tools.hpp>
#include <gsl/gsl_rng.h>
gsl_rng* global_rng;

const double FMGEV = 5.08;

const int INTPOINTS = 3;
const double INTACCURACY = 0.01;

double r0 = 0.3 * FMGEV;    // Dipole size used in the analysis
double N0 = 0.03;

struct protonsizehelper
{
    IPGlasma *glasma;
    double b;
    double theta_b;
    double r;
};

double AngularAverageAmplitude_theta_r(double theta_r, void* p)
{
    protonsizehelper* par = (protonsizehelper*) p;
    double x1[2] = { par->b * cos(par->theta_b) + 0.5*par->r * cos(theta_r), par->b * sin(par->theta_b) + 0.5*par->r * sin(theta_r) };
    double x2[2] = { par->b * cos(par->theta_b) - 0.5*par->r * cos(theta_r), par->b * sin(par->theta_b) - 0.5*par->r * sin(theta_r) };
    
    return par->glasma->Amplitude(0.01, x1, x2);
    
}

// Average over dipole orientation
double AverageAmplitude_theta_b(double theta_b, void* p)
{
    protonsizehelper* par = (protonsizehelper*) p;
    par->theta_b=theta_b;
    
    // Do angular integral
    gsl_function fun;
    fun.function =AngularAverageAmplitude_theta_r;
    fun.params = par;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, 1e-5, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    
    return result / (2.0*M_PI);
}

// Average over impact parameter angle
double AverageAmplitude_b(double b, void* p)
{
    protonsizehelper* par = (protonsizehelper*) p;
    par->b=b;
    
    // Do angular integral
    gsl_function fun;
    fun.function =AverageAmplitude_theta_b;
    fun.params = par;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, 1e-5, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    
    //cout << "b=" << b << " avg " <<result / (2.0 * M_PI) << endl;
    
    return result / (2.0 * M_PI)- N0;

    
}




int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename
    string fname = argv[1];
    N0 = StrToReal(argv[2]);    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    IPGlasma glasma(fname);
    protonsizehelper par;
    par.r = r0;
    par.glasma = &glasma;
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    double x_lo=0.001;
    double x_hi=10;
    gsl_function F;
    F.function = &AverageAmplitude_b;
    F.params = &par;
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    double r;
	
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.001);
        
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    cout << r << endl;
    
/*    
    for (double b=0; b<15; b+=0.1)
    {
        double average_n = AverageAmplitude_b(b, &par);
        cout << b << " " << average_n << endl;
    }
    
 */   
    
    
}
