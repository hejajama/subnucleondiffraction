/* Calculate avarage dipole amplitude averaged over b
*/

#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <tools/tools.hpp>
#include <gsl/gsl_rng.h>

gsl_rng* global_rng;
using namespace std;

struct bhelper
{
    IPGlasma *glasma;
    double r;
    double theta_r;
    double theta_b;
    double b;
};

const int INTPOINTS = 10000;
const double INTACCURACY = 0.05;
double rhelperf(double theta_r, void* p); // average over dipole orientation
double bhelperf(double theta_b, void* p); // average over b angle

int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b
    string fname = argv[1];
    double b= StrToReal(argv[2]);
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    IPGlasma glasma(fname);
    bhelper helper;
    helper.glasma = &glasma;
    helper.b=b;
    
    Ipsat_Proton ipsat;
    ipsat.SetProtonWidth(0);
    ipsat.SetQuarkWidth(4.0);
    ipsat.InitializeTarget();
    
    gsl_function f;
    f.function = &bhelperf;  // do a circle around the proton
    f.params = &helper;
    
    for (double r=1e-2; r<10; r*=1.2)
    {
        helper.r=r;
    
         
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
        double result,error;
        int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        if (status)
            cerr << "#btheta_b int failed, result  " << result << " relerror " << error << " r " <<r  << endl;
        
        
        gsl_integration_workspace_free(w);
        
        result = result / (2.0*M_PI*2.0*M_PI);
        
        // IPsat comparison
        Vec q1(b+0.5*r,0); Vec q2(b-0.5*r,0);
        double ipsat_n = ipsat.DipoleAmplitude::Amplitude(1e-3, q1, q2);
        
        cout << r << " " << result << " " << ipsat_n << endl;
    }
    
    


}


double bhelperf(double theta_b, void* p)
{
    bhelper* par = (bhelper*)p;
    par->theta_b = theta_b;
    gsl_function f;
    f.function = &rhelperf;  // do a circle around the proton
    f.params = par;
    
    
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#btheta_r int failed, result  " << result << " relerror " << error << " theta_b " <<theta_b << endl;
    
    gsl_integration_workspace_free(w);
    
    return result;

    
}

double rhelperf(double theta_r, void* p)
{
    bhelper* par = (bhelper*)p;
    par->theta_r = theta_r;
    
    Vec b(par->b * cos(par->theta_b), par->b * sin(par->theta_b));
    Vec r(par->r * cos(par->theta_r), par->r * sin(par->theta_r));
    Vec r2 = r*0.5;
    Vec q1 = b + r2;
    Vec q2 = b - r2;
    
    return par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
}