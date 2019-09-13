/* Calculate avarage dipole amplitude averaged over b
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
    double theta_r;
    double theta_b;
    double b;
};

const int INTPOINTS = 2;
const double INTACCURACY = 0.5;
const double MINB=1;
double MAXB=1.25;
double rhelperf(double theta_r, void* p); // average over dipole orientation
double bhelperf(double b, void* p); // average over b
double bhelperf_theta(double theta_b, void* p);

int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    double rc= StrToReal(argv[3]);
    double b = StrToReal(argv[2]);
    double step = StrToReal(argv[4]);
    cout << "# Filename: " << fname << " b " << b << " Gev^-1 schwinger rc " << rc << " step[fm]" << endl;

    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    IPGlasma glasma(fname, step);
	if (rc > 0)
	glasma.SetSchwinger(true, rc);
    bhelper helper;
    helper.glasma = &glasma;
	helper.b=b; helper.theta_b=0;
    //helper.b=b;
    
    Ipsat_Proton ipsat;
    ipsat.SetProtonWidth(0);
    ipsat.SetQuarkWidth(4.0);
    ipsat.InitializeTarget();
    
    gsl_function f;
    f.function = &bhelperf;  // do a circle around the proton
    //f.function = &rhelperf;
	f.params = &helper;
    
    cout << "# r [1/GeV]  N(x=" << argv[2] << ", y=+/- r/2),  <N(r, " << b+MINB << "<b<"<<b+MAXB << "), IPsat N(r,b=" << argv[2] << ")" << endl;
    
    for (double r=5e-3; r<40; r*=1.1)
    {
        double x1[2]={StrToReal(argv[2]),r/2.0}; double x2[2]={StrToReal(argv[2]),-r/2.0};
        double fixedresult = glasma.Amplitude(1e-3,x1,x2);
    
        helper.r = r;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
	//MAXB = (5.12/2.0*5.068-r/2.0);	// so that q/qbar is never outside the lattice
        double result,error;
        int status = gsl_integration_qag(&f, b+MINB, b+MAXB, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
	//int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 1e-5, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        //if (status)
          //  cerr << "#bint failed (" << gsl_strerror (status) << ", result  " << result << " relerror " << error/result << " r " <<r  << endl;
        
        
        gsl_integration_workspace_free(w);
  		//result /= (2.0*M_PI);     
        result = result / (2.0*M_PI*M_PI*( (MAXB + b)*(MAXB+b) - (MINB+b)*(MINB+b)) );
	//	result = result / (2.0*M_PI*2.0*M_PI);
        
        // IPsat comparison
	double b = StrToReal(argv[2]);
        Vec q1(b+0.5*r,0); Vec q2(b-0.5*r,0);
        double ipsat_n = ipsat.DipoleAmplitude::Amplitude(1e-2, q1, q2);
        
        cout << r << " " << fixedresult << " " << result << " " << ipsat_n << endl;
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
    
   	// Fix r angle by theta_b
   	// r parallel to b
	double theta_r = theta_b; //+M_PI/2.0;
   	return rhelperf(theta_r, par); 

    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    
   // if (status)
     //   cerr << "#btheta_r int failed, result  " << result << " relerror " << error << " theta_b " <<theta_b << endl;
    
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
    //Vec q1 = b; // + r2;
    //Vec q2 = b + r; //- r2;
    Vec q1 = b + r2;
    Vec q2 = b - r2;
    
    return par->glasma->DipoleAmplitude::Amplitude(0.01, q1, q2);
}
