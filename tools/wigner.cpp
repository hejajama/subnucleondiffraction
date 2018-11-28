/*
* Calculate avarage dipole amplitude averaged over b
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
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>


gsl_rng* global_rng;
using namespace std;

struct inthelper_wigner
{
    DipoleAmplitude* dipole;
    double theta_b;
    double b;
    double q;
    bool real_part;
    double xpom;
    bool ipglasma;  // if true, FT N, if False, FT S (dies at large r)
};

double inthelperf_mc( double *vec, size_t dim, void* par);
int MCINTPOINTS_WIGNER = 2e7;
int main(int argc, char* argv[])
{
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();

    inthelper_wigner helper;
    helper.real_part = true;
    double b=2.5;
    helper.xpom=0.01;
    
    double A=1;
    
    std::string ipglasmafile="";
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-b")
            b=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-q")
            helper.q = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-mcintpoints")
            MCINTPOINTS_WIGNER = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-A")
            A=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xpom")
            helper.xpom=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-imag")
            helper.real_part=false;
        else if (string(argv[i])=="-ipglasma")
            ipglasmafile = argv[i+1];
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
    }
                
    helper.b=b;
    
    DipoleAmplitude* dipole;
    if (ipglasmafile=="")
    {
        dipole = new Ipsat_Proton;
        ((Ipsat_Proton*)dipole)->SetProtonWidth(0);
        ((Ipsat_Proton*)dipole)->SetQuarkWidth(4);
        ((Ipsat_Proton*)dipole)->SetA(1);
        ((Ipsat_Proton*)dipole)->InitializeTarget();
        helper.ipglasma=false;
    }
    else
    {
        dipole = new IPGlasma(ipglasmafile, 0.00731429, BINARY);
        helper.ipglasma=true;
    }
    //IPGlasma dipole(argv[1], 0.00731429, BINARY);

    helper.dipole = dipole;
    
    double *lower, *upper;
    
    lower = new double[3];
    upper = new double[3];
    
    
    
    lower[0]=lower[1]=0;
    upper[0] = 10*5.068 ; // Max r
    upper[1] = 2.0*M_PI;
    
    
    gsl_monte_function F;
    F.f = &inthelperf_mc;
    
    F.dim=2;
    if (helper.ipglasma==true)
    {
        F.dim = 3;
        lower[2]=0;
        upper[2]=2.0*M_PI;   // Overall rotation
    }
        
    
    F.params = &helper;
    
    double result,error;
    
     gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        cout << "# b = " << helper.b << " q = " << helper.q << " A = " << A << " xp = " << helper.xpom  <<  endl;
    cout << "# angle  wigner  montecarloerror" << endl;
   // for (double k=0.4; k<15; k*=1.5){
        
        for (double th = 0; th<= 2.0*M_PI*1.0001; th += 2.0*M_PI/30)
        {
            helper.theta_b = th;
            
            GRADIENT = false;
            double r = 0;
            double laplace=0;
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &r, &error);
            
            if (helper.ipglasma == true) // IPGlasma only currenltly repsects global GRADIENT
            {
                GRADIENT = true;
                GRADIENT_STEP=1;
                
                gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &laplace, &error);
            }
            else
            {
            
            
                 // Laplace operator:
                 // nabla_b^2 = \partial_r^2 + 1/r \partial_r + 1/r^2 \partial_th
                 // First try to use simple finite difference method
                 // d^2 f/dr^2 = [f(x+h) + f(x-h) -f(x)]/h^2
                 double db = 0.2*5.068; // Seems to make sense, quite large, but hopefully
                 // laplace term is a small correction anyway
                
                
                
                 if (b - db < 0)
                 {
                 cerr<< "Computing laplace only works for b>0" << endl;
                 exit(1);
                 }
                double b_plus;
                helper.b = b + db;
                gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &b_plus, &error);
                
                double b_minus;
                helper.b = b - db;
                gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &b_minus, &error);
                
                 double delta_th = 0.15;
                 if (th - delta_th < 0 or th + delta_th > 2.0*M_PI)
                 {
                     cerr << "Computing laplace only works if angle is not 0 or 2pi" << endl;
                     exit(1);
                 }
                helper.theta_b = th + delta_th;
                double th_plus;
                gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &th_plus, &error);
                helper.theta_b = th - delta_th;
                double th_minus;
                gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_WIGNER, global_rng, s, &th_minus, &error);
                
                 double der2b = (b_plus + b_minus - 2.0*r)/(db*db);
                 double derb = (b_plus - b_minus)/(2.0*db);
                 double der2th = (th_plus + th_minus -2.0*r)/(delta_th*delta_th);
                
                 laplace = der2b + 1.0/b * derb  + 1.0/(b*b)*der2th;
            }
            
            if (helper.ipglasma==true)
            {
                // Azimuthal average
                r /= 2.0*M_PI;
                laplace /= 2.0*M_PI;
            }
            cout << th << " " << r << " " << laplace << endl;
            
            

        }
    
   
    
    gsl_monte_miser_free(s);
   
    


}


double inthelperf_mc( double *vec, size_t dim, void* p)
{
    inthelper_wigner* par = (inthelper_wigner*)p;
    double r = vec[0];
    double theta_r = vec[1];
   
    
    double bx = par->b*cos(par->theta_b);
    double by = par->b*sin(par->theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    
    double b = par->b;
    double theta_b = par->theta_b;
   
    DipoleAmplitude* dipole = par->dipole;
    
    double z=0.5;
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    
    double q1[2] = {qx,qy};
    double q2[2] = {qbarx, qbary};
    

    double b_dot_q = bx*rx+by*ry;
    
    std::complex<double> imag(0,1);
    complex<double> exponent = imag*b_dot_q;
    
    // Overall rotation
    if (par->ipglasma)
    {
        double overall_angle = vec[2];
        q1[0] = std::cos(overall_angle)*q1[0] + std::sin(overall_angle)*q1[1];
        q1[1] = -std::sin(overall_angle)*q1[0] + std::cos(overall_angle)*q1[1];
        
        q2[0] = std::cos(overall_angle)*q2[0] + std::sin(overall_angle)*q2[1];
        q2[1] = -std::sin(overall_angle)*q2[0] + std::cos(overall_angle)*q2[1];
    }
    double amp = dipole->Amplitude(par->xpom,q1,q2);

    
    
    if (par->ipglasma==false)   // IPsat, FT S; not N
        amp = 1.0 - amp;
    
    complex<double> result = std::exp(exponent) * amp;
    
    // Prefactors and Jacobian
    result *= r;
    // Todo nc, as, 2pi
    
    if (par->real_part)
        return result.real();
    else
        return result.imag();
    return 0;
}
