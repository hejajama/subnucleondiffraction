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

struct inthelper_husimi
{
    DipoleAmplitude* dipole;
    double theta_b;
    double b;
    double k;
    double l;
    bool real_part;
    double xpom;
   
};

double inthelperf_mc( double *vec, size_t dim, void* par);
int MCINTPOINTS_HUSIMI = 2e7;
bool avereages_azimuth = false;

int main(int argc, char* argv[])
{
    // Params ipglasmafname k theta
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();

    inthelper_husimi helper;
    helper.real_part = true;
    helper.k = StrToReal(argv[2]);
    helper.b=2.5;
    helper.l=1;
    helper.xpom=0.01;
    
    double A=1;
    
    std::string ipglasmafile="";
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-b")
            helper.b=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-k")
            helper.k = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-l")
            helper.l = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-mcintpoints")
            MCINTPOINTS_HUSIMI = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-A")
            A=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xpom")
            helper.xpom=StrToReal(argv[i+1]);
        else if (string(argv[i])=="-imag")
            helper.real_part=false;
        else if (string(argv[i])=="-ipglasma")
            ipglasmafile = argv[i+1];
        else if (string(argv[i])=="-azimuth_averages")
            avereages_azimuth =true;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
    }
                
    
    
    DipoleAmplitude* dipole;
    if (ipglasmafile=="")
    {
        dipole = new Ipsat_Proton;
        ((Ipsat_Proton*)dipole)->SetProtonWidth(0);
        ((Ipsat_Proton*)dipole)->SetQuarkWidth(4);
        ((Ipsat_Proton*)dipole)->SetA(1);
        ((Ipsat_Proton*)dipole)->InitializeTarget();
    }
    else
    {
        dipole = new IPGlasma(ipglasmafile, 0.00731429, BINARY);
    }
    //IPGlasma dipole(argv[1], 0.00731429, BINARY);

    helper.dipole = dipole;
    
    double *lower, *upper;
    
    lower = new double[4];
    upper = new double[4];
    
    
    
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 10*5.068 ; // Max r
    upper[1] = 2.0*M_PI;
    upper[2] = 10*5.068;    // max b
    upper[3] = 2.0*M_PI;
    
    gsl_monte_function F;
    F.dim = 4;
    
    if (avereages_azimuth == true)
    {
        F.dim=5;    // Additional dimension: overall rotation
        lower[4]=0;
        upper[4]=2.0*M_PI;
    }
    
    F.f = &inthelperf_mc;
    
    
    F.params = &helper;
    
    double result,error;
    
     gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        cout << "# b = " << helper.b << " k = " << helper.k << " A = " << A << " xp = " << helper.xpom  << " l= " << helper.l << endl;
    cout << "# angle  Husimi  montecarloerror" << endl;
   // for (double k=0.4; k<15; k*=1.5){
        
        for (double th = 0; th<= 2.0*M_PI*1.0001; th += 2.0*M_PI/30)
        {
           
            
            helper.theta_b = th;
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_HUSIMI, global_rng, s, &result, &error);
                //cout <<"# rot " << helper.overall_angle << " res " << result <<  " pm " << error <<endl;
            
            cout << th << " " << result/(2.0*M_PI) << " " << error/(2.0*M_PI) << endl;
        }
        /*std::vector<double> res;
        
        for (double th = 0; th <= M_PI/2.0*1.01; th += M_PI/2.0)
        {
            helper.k = k;
            helper.theta_b = th;
            helper.real_part=true;
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
            cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
            res.push_back(result);
            
            helper.real_part=false;
            result=0;
            //gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
            //cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
            res.push_back(result);
            
        }
        cout << k << " " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << endl;
         */
        
        
 //   }
    

   
    
    gsl_monte_miser_free(s);
   
    


}


double inthelperf_mc( double *vec, size_t dim, void* p)
{
     inthelper_husimi *par = (inthelper_husimi*)p;
    
    double r = vec[0];
    double theta_r = vec[1];
    double b2 = vec[2]; // b'
    double theta_b2 = vec[3];
    double theta_b = par->theta_b;
    

    
    
    double b2x = b2*cos(theta_b2);
    double b2y = b2*sin(theta_b2);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    
   
    double b = par->b;
    
    double k = par->k;
    double l = par->l;
    DipoleAmplitude* dipole = par->dipole;
    
    double z=0.5;
    double qx = b2x + z*rx; double qy = b2y + z*ry;
    double qbarx = b2x - (1.0-z)*rx; double qbary = b2y - (1.0-z)*ry;
    
    double q1[2] = {qx,qy};
    double q2[2] = {qbarx, qbary};
    
    
    double b_minus_b2_sqr = b*b + b2*b2 - 2.0*b*b2*std::cos(theta_b - theta_b2);
    
    // (k + i/(2l)*r)^2
    std::complex<double> imag(0,1);
    std::complex<double> kilr_sqr = k*k - std::pow(r/(2*l), 2.0) + imag/l * k*r*std::cos(theta_r);
    
    complex<double> exponent = -1.0/(l*l) * b_minus_b2_sqr - r*r/(4.0*l*l) + imag*k*r*std::cos(theta_r);
    
    if (avereages_azimuth == true)
    {// Rotate dipole
        double overall_angle = vec[4];
        q1[0] = std::cos(overall_angle)*q1[0] + std::sin(overall_angle)*q1[1];
        q1[1] = -std::sin(overall_angle)*q1[0] + std::cos(overall_angle)*q1[1];
        
        q2[0] = std::cos(overall_angle)*q2[0] + std::sin(overall_angle)*q2[1];
        q2[1] = -std::sin(overall_angle)*q2[0] + std::cos(overall_angle)*q2[1];
    }
    complex<double> result = std::exp(exponent) * (1.0/(l*l)*b_minus_b2_sqr + l*l*kilr_sqr)
        * dipole->Amplitude(par->xpom,q1,q2);
    
    // Prefactors and Jacobian
    result *= -b2*r/(l*l);
    // Todo nc, as, 2pi
    
    if (par->real_part)
        return result.real();
    else
        return result.imag();
    return 0;
}
