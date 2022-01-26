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
    double theta_r;
    double theta_b;
    double b;
};

const int INTPOINTS = 20;
const double INTACCURACY = 0.01;
double rhelperf(double theta_r, void* p); // average over dipole orientation
double bhelperf(double b, void* p); // average over b
double bhelperf_theta(double theta_b, void* p);

int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];

    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    IPGlasma glasma(fname, 0.01, BINARY );
    bhelper helper;
    helper.glasma = &glasma;
	helper.theta_b=0;
    //helper.b=b;
    
       gsl_function f;
    f.function = &bhelperf_theta; // Rotate r
    //f.function = &rhelperf;
	f.params = &helper;
   
    double r = 0.2 * FMGEV;  
    double MAXB = 4.5 * FMGEV;;
    cout <<"# Dipole size r=" << r/FMGEV << " fm" << endl;
    cout <<"# Datafile: " << fname << endl;
    for (double b = -MAXB; b < MAXB; b+= 0.05) 
    {
        helper.b=b;
        double x1[2]={b-r/2.0,0}; double x2[2]={b+r/2.0,0};
        double fixedresult = glasma.Amplitude(1e-3,x1,x2);
    
        helper.r = r;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
        double result,error;
	    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 1e-5, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        //if (status)
          //  cerr << "#bint failed (" << gsl_strerror (status) << ", result  " << result << " relerror " << error/result << " r " <<r  << endl;
        
        
        gsl_integration_workspace_free(w);
  		result /= (2.0*M_PI);     
        
        
        cout << b << " " << fixedresult << " " << result <<  endl;
    }
    
    


}


double bhelperf_theta(double theta_r, void* p)
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
