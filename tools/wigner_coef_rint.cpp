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

const int INTPOINTS = 2;
const double INTACCURACY = 0.01;
double helperf_theta(double theta_b, void* p);
double helperf_r(double r, void* p);

int wigner_coef = 0;     //1=xW1, 0 = xW0

gsl_integration_workspace *w_overall;
gsl_integration_workspace *w_relative;

int main(int argc, char* argv[])
{
    const double fmgev=5.068;
    // Arguments: ipglasma filename step  cos2phimom
    if (argc != 3)
    {
        cerr << "Arguments: filaneme wlinestepsize [fm] " << endl;
        return 1;
    }
    string fname = argv[1];
    double step = StrToReal(argv[2]);
    cout << "# Filename: " << fname <<  "  step[fm] " << step <<  endl;
    double step_gev = step *fmgev;
    double maxb = 1.11*fmgev;
    double bstep = 0.05 * fmgev;
    double maxk = 4.01;
    double kstep = 0.25;
    
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
    
    cout << "#k  b   integral" << endl;
    for (double k=0; k <=maxk;  k+=kstep)
    {
        for (double b=0; b <= maxb; b+=bstep)
        {
            
            helper.k=k;
            helper.b=b;
            
            gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
            double result,error;
            int status = gsl_integration_qag(&f, 0, 4*fmgev, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
            
            
            gsl_integration_workspace_free(w);
            
            
            result /= (2.0*M_PI);
            cout << k << " " << b << " " << result << endl;
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


