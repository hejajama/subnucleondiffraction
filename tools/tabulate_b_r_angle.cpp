/* Calculate avarage dipole amplitude averaged over overall azimuthal angle of b, keeping angle between b and r fixed 
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

const int INTPOINTS = 3;
const double INTACCURACY = 0.1;
double bhelperf_theta(double theta_b, void* p);

int main(int argc, char* argv[])
{
	const double fmgev=5.068;
    // Arguments: ipglasma filename step 
    if (argc != 3)
    {
        cerr << "Arguments: filaneme wlinestepsize [fm] " << endl;
        return 1;
    }
    string fname = argv[1];
    double step = StrToReal(argv[2]);
    cout << "# Filename: " << fname <<  "  step[fm] " << step <<  endl; 
    double step_gev = step *fmgev;
    double maxr = 3.0*fmgev;
    double maxb = 3.0*fmgev;
    int points_r = 300;
    int points_b = 300;

    if (maxr / points_r < step_gev)
    {
         cerr << "Step size is less than grid spacing!" << endl;
         exit(1);
    }
    
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
    f.function = &bhelperf_theta;  // do a circle around the proton
	f.params = &helper;
    
   cout << "#r  b   theta_r_b   N" << endl; 
    for (double r=0; r <=maxr;  r+=maxr/points_r)
    {
	for (double b=0; b <= maxb; b+=maxb / points_b)
	{
    		for (double theta_r_b=0; theta_r_b <= M_PI/2.0; theta_r_b += (M_PI/2.0)/20.0)
   		 {
    
        		helper.theta_r_b = theta_r_b;
			helper.r=r;
			helper.b=b;
        		gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
        		double result,error;
			int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 1e-7, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        
       			 gsl_integration_workspace_free(w);
  			result /= (2.0*M_PI);     
        		cout << r << " " << b << " " << theta_r_b << " " << result << endl;
		}
	}    
    }
    
    
	return 0;

}



double bhelperf_theta(double theta_b, void* p)
{
    bhelper* par = (bhelper*)p;

    Vec b(par->b * cos(theta_b), par->b * sin(theta_b));
    Vec r(par->r * cos(theta_b+par->theta_r_b), par->r * sin(theta_b + par->theta_r_b));
    
   Vec r2 = r*0.5;
    //Vec q1 = b; // + r2;
    //Vec q2 = b + r; //- r2;
    Vec q1 = b + r2;
    Vec q2 = b - r2;
   
    return par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
 

    
}


