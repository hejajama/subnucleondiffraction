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
    bool imaginary_part;
};

const int INTPOINTS = 4;
const double INTACCURACY = 0.015;
double MINB=0;
double MAXB=3*5.068;
double rhelperf(double theta_r, void* p); // average over dipole orientation
double bhelperf(double b, void* p); // average over b
double bhelperf_theta(double theta_b, void* p);

int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    cout << "# Filename: " << fname << endl;

    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    IPGlasma glasma(fname, 0.01, BINARY);
    bhelper helper;
    helper.glasma = &glasma;
   glasma.SetPeriodicBoundaryConditions(true); 
    gsl_function f;
    f.function = &bhelperf;  // do a circle around the proton
    //f.function = &rhelperf;
	f.params = &helper;
    
    double fmgev=5.068;
    
    for (double r=0.01*fmgev; r<2*fmgev; r+=0.01*fmgev)
    {
       // MINB=0; MAXB=2.0*M_PI;
        helper.r = r;
        helper.b=0; helper.theta_b=0;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
        helper.imaginary_part=false;
        double result,error;
        int status = gsl_integration_qag(&f, MINB, MAXB, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
        

        helper.imaginary_part=true;
        double result_i,error_i;
        int status_i = gsl_integration_qag(&f, MINB, MAXB, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result_i, &error_i);
 
        //if (status)
          //  cerr << "#bint failed (" << gsl_strerror (status) << ", result  " << result << " relerror " << error/result << " r " <<r  << endl;
        
        
        gsl_integration_workspace_free(w);
  //		result /= (2.0*M_PI);    
  //      result_i /= (2.0*M_PI); 
        result = result / (M_PI*( (MAXB )*(MAXB) - (MINB)*(MINB)) );
	result_i = result_i / (M_PI*( (MAXB )*(MAXB) - (MINB)*(MINB)) );

        // FIxed configuration for comparison
        
        double q1[2] = {-r, -r};
        double q2[2] = {r,-r};
        double q3[2] = {0, (std::sqrt(3)-1.0)*r};
        
        std::complex<double> fixed = glasma.BaryonOperator(0.01, q1,q2,q3);

	//	result = result / (2.0*M_PI*2.0*M_PI);
        
        
        cout << r << " " << result << " " << result_i << " " << fixed.real() << " " << fixed.imag() << endl;
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
    double r = par->r;


/*
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, INTACCURACY, INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &error);
    
   // if (status)
     //   cerr << "#btheta_r int failed, result  " << result << " relerror " << error << " theta_b " <<theta_b << endl;
    
    gsl_integration_workspace_free(w);
    
    return result;
*/



    Vec b(par->b*std::cos(theta_b), par->b*std::sin(theta_b));
    Vec q1v (-par->r, -par->r);
    Vec q2v (par->r, -par->r);
    Vec q3v (0, (std::sqrt(3)-1.0)*par->r);
    q1v = b + q1v;
    q2v = b + q2v;
    q3v = b + q3v;

    double q1[2] = {q1v.GetX(), q1v.GetY()} ;
    double q2[2] = {q2v.GetX(), q2v.GetY()} ;
    double q3[2] = {q3v.GetX(), q3v.GetY()} ;

    std::complex<double> res = par->glasma->BaryonOperator(0.01, q1, q2, q3); 
    if (par->imaginary_part) return res.imag(); 
    else return res.real();
    
}

double rhelperf(double theta_r, void* p)
{
    bhelper* par = (bhelper*)p;
    par->theta_r = theta_r;


    Vec b(par->b*std::cos(par->theta_b), par->b*std::sin(par->theta_b));
    Vec q1v (-par->r, -par->r);
    Vec q2v (par->r, -par->r);
    Vec q3v (0, (std::sqrt(3)-1.0)*par->r);
    q1v.Rotate2D(theta_r);
    q2v.Rotate2D(theta_r);
    q3v.Rotate2D(theta_r);
    q1v = b + q1v;
    q2v = b + q2v;
    q3v = b + q3v;

    double q1[2] = {q1v.GetX(), q1v.GetY()} ;
    double q2[2] = {q2v.GetX(), q2v.GetY()} ;
    double q3[2] = {q3v.GetX(), q3v.GetY()} ;

    std::complex<double> res = par->glasma->BaryonOperator(0.01, q1, q2, q3); 
    if (par->imaginary_part) return res.imag(); 
    else return res.real();
 

}

