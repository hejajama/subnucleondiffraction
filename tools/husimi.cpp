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
#include <gsl/gsl_sf_bessel.h>


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
    
    double r;
    gsl_integration_workspace *ws1;
    gsl_integration_workspace *ws2;
   
};

const int INTEGRATION_INTERVALS = 8;
const double COS2PHIACCURACY = 0.001;
double inthelperf_mc( double *vec, size_t dim, void* par);
double inthelperf_mc_xH1( double *vec, size_t dim, void* par);
int MCINTPOINTS_HUSIMI = 2e7;
bool avereages_azimuth = false;
bool COMPUTE_HUSIMI_V2=false;   // Compute xH1 as in 1609.05773 eq (26)
bool V0=false;


double inthelperf_cos2phi_n(double phi_rb, void* p)
{
    inthelper_husimi *par = (inthelper_husimi*)p;
    
    double rx = par->r * cos(par->theta_b + phi_rb);
    double ry = par->r * sin(par->theta_b + phi_rb);
    
    double bx = par->b * cos(par->theta_b);
    double by = par->b * sin(par->theta_b);
    
    double z=0.5;
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    
    double q1[2] = {qx,qy};
    double q2[2] = {qbarx, qbary};
    
    return std::cos(2.0*phi_rb) * par->dipole->Amplitude(0.01, q1,q2);

}


double inthelperf_cos2phi_overall(double overall_rotation, void* p)
{
    inthelper_husimi *par = (inthelper_husimi*)p;
    par->theta_b = overall_rotation;
    
    double result,abserr;
    gsl_function f;
    f.function = &inthelperf_cos2phi_n;
    f.params = par;
    
    int status = gsl_integration_qag(&f, 0, 2*M_PI, 0, COS2PHIACCURACY, INTEGRATION_INTERVALS, GSL_INTEG_GAUSS51, par->ws2, &result, &abserr);
    
    return result/(2.0*M_PI);
    
}




double Cos2Phi_N(double r, double b, DipoleAmplitude *dipole)
{
    // Compute \int_0^2pi cos(2phi) N(r,b,phi)
    
    inthelper_husimi par;
    par.b=b;
    par.r=r;
    par.dipole = dipole;
    par.ws1 = gsl_integration_workspace_alloc(INTEGRATION_INTERVALS);
    par.ws2 = gsl_integration_workspace_alloc(INTEGRATION_INTERVALS);
    
    double result,abserr;
    gsl_function f;
    f.function = &inthelperf_cos2phi_overall;
    f.params = &par;
    
    
    int status = gsl_integration_qag(&f, 0, 2*M_PI, 0, COS2PHIACCURACY, INTEGRATION_INTERVALS, GSL_INTEG_GAUSS51, par.ws1, &result, &abserr);
    
    if (status)
    {
        //cerr << "Cos2phi integral failed (r=" << r <<", b=" << b <<"), result " << result << " +/- " << abserr << endl;
    }
    //cout << r << " " << b << " " << result << " " << abserr << endl;
    
    
    delete par.ws1;
    delete par.ws2;
    
    return result;
    
}

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
        else if (string(argv[i])=="-xH1")
            COMPUTE_HUSIMI_V2 = true;
        else if (string(argv[i])=="-v0")
            V0 = true;
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
        dipole = new IPGlasma(ipglasmafile, 0.005, BINARY); //0.00731429, BINARY);
    }
    //IPGlasma dipole(argv[1], 0.00731429, BINARY);

    helper.dipole = dipole;
    
    double *lower, *upper;
    
    lower = new double[5];
    upper = new double[5];
    
    gsl_monte_function F;
    if (COMPUTE_HUSIMI_V2 == false)
    {
        lower[0]=lower[1]=lower[2]=lower[3]=0;
        upper[0] = 8*5.068 ; // Max r
        upper[1] = 2.0*M_PI;
        upper[2] = 10*5.068;    // max b
        upper[3] = 2.0*M_PI;
        
        
        F.dim = 4;
        
        if (avereages_azimuth == true)
        {
            F.dim=5;    // Additional dimension: overall rotation
            lower[4]=0;
            upper[4]=2.0*M_PI;
        }
        
        F.f = &inthelperf_mc;
    }
    else
    {
        
        lower[0]=lower[1]=lower[2]=lower[3]=0;
        upper[0]=5*5.068; // max r
        upper[1]=5*5.068; // max b2
        upper[2] = 2.0*M_PI;   // r,b2 angle
        upper[3] = 2.0*M_PI;    // overall rotation
        
        F.dim=4;
        
        /*
        lower[0]=std::log(1e-3);
        lower[1]=0;
        upper[0]=std::log(50);
        upper[1]=3*5.068; // max b2
        F.dim=2;
         */
        
        F.f = &inthelperf_mc_xH1;
        
        avereages_azimuth = true;  // This is done automatically in this case, true = correct normalization?
        
    }
    
    
    F.params = &helper;
    
    double result,error;
    
     gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        cout << "# b = " << helper.b << " k = " << helper.k << " A = " << A << " xp = " << helper.xpom  << " l= " << helper.l << endl;
   
    
    if (COMPUTE_HUSIMI_V2 == false)
    {

         cout << "# angle  Husimi  montecarloerror" << endl;
        for (double th = 0; th<= 2.0*M_PI*1.0001; th += 2.0*M_PI/30)
        {
           
            
            helper.theta_b = 0; //th;
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_HUSIMI, global_rng, s, &result, &error);
                //cout <<"# rot " << helper.overall_angle << " res " << result <<  " pm " << error <<endl;
            
            if (avereages_azimuth)
            {
                result /= 2.0*M_PI;
                error /= 2.0*M_PI;
            }
            
            
            cout << th << " " << result << " " << error << endl;
        }
    }
    else
    {
         cout << "# k  Husimi  montecarloerror" << endl;
	double kstep = 0.25;
        //for (double k=0.5; k<3.51; k += kstep) //3.6; k+=0.2)
	for (double k=0.1; k<0.51; k+=0.1)
        {
            helper.k=k;
            
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS_HUSIMI, global_rng, s, &result, &error);
            result /= 2.0*M_PI;
            error /= 2.0*M_PI;
            
            cout << k << " " << result << " " << error << endl;

//            if (k>0.15) kstep = 0.2;
        }
        
    }
        
        
 //   }
    

        
           // for (double k=0.4; k<15; k*=1.5){
   
    
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
    
    // (k + i/(2l^2)*r)^2
    std::complex<double> imag(0,1);
    std::complex<double> kilr_sqr = k*k - std::pow(r/(2*l*l), 2.0) + imag/(l*l) * k*r*std::cos(theta_r);
    
    complex<double> exponent = -1.0/(l*l) * b_minus_b2_sqr - r*r/(4.0*l*l) + imag*k*r*std::cos(theta_r);
    
    double overall_angle = vec[4];
    double rot_q1[2];
    double rot_q2[2];
    
    double amplitude=0;
    
    if (avereages_azimuth == true)
    {// Rotate dipole
        cout << "hi" << endl;
        rot_q1[0] = std::cos(overall_angle)*q1[0] + std::sin(overall_angle)*q1[1];
        rot_q1[1] = -std::sin(overall_angle)*q1[0] + std::cos(overall_angle)*q1[1];
        
        rot_q2[0] = std::cos(overall_angle)*q2[0] + std::sin(overall_angle)*q2[1];
        rot_q2[1] = -std::sin(overall_angle)*q2[0] + std::cos(overall_angle)*q2[1];
        
        amplitude =dipole->Amplitude(par->xpom,rot_q1,rot_q2);
    }
    else
    {
        amplitude =dipole->Amplitude(par->xpom,q1,q2);
    }
    complex<double> result = std::exp(exponent) * (1.0/(l*l)*b_minus_b2_sqr + l*l*kilr_sqr) * amplitude;
    
    // Prefactors and Jacobian
    result *= -b2*r/(l*l);
    // Todo nc, as, 2pi
    
    if (par->real_part)
        return result.real();
    else
        return result.imag();
    return 0;
}




double inthelperf_mc_xH1( double *vec, size_t dim, void* p)
{
    inthelper_husimi *par = (inthelper_husimi*)p;
    
    /*
    double r = std::exp(vec[0]);
    double b2 = vec[1]; // b'
     */
    
    double r = vec[0];
    double b2 = vec[1];
    double theta_rb2 = vec[2];
    double overall_rotation = vec[3];
    
    double b2x = b2*cos(overall_rotation);
    double b2y = b2*sin(overall_rotation);
    
    double rx = r*cos(overall_rotation + theta_rb2);
    double ry = r*sin(overall_rotation + theta_rb2);
    
    
  
    
    double b = par->b;
    
    double k = par->k;
    double l = par->l;
    DipoleAmplitude* dipole = par->dipole;
    

    
    double z=0.5;
    double qx = b2x + z*rx; double qy = b2y + z*ry;
    double qbarx = b2x - (1.0-z)*rx; double qbary = b2y - (1.0-z)*ry;
    
    double q1[2] = {qx,qy};
    double q2[2] = {qbarx, qbary};
    
    
    double result = 0;
    
    double i2 = gsl_sf_bessel_In(2, 2.0*b*b2/(l*l));
    double i1 = gsl_sf_bessel_I1(2*b*b2/(l*l));
    double j2 = gsl_sf_bessel_Jn(2,k*r);
    double j1 = gsl_sf_bessel_J1(k*r);
    
    
    if (V0 == false)
    {
        result = 1.0/(2.0*M_PI) * std::exp(-1.0/(l*l) *(b*b+b2*b2) - r*r/(4.0*l*l));
        
        result *=(((1.0/(l*l) * (b*b + b2*b2) + l*l*k*k - r*r/(4.0*l*l)) * i2 - 2.0*b*b2 /(l*l) * i1 ) * j2  + k*r*i2*j1)
            * std::cos(2.0*theta_rb2) * dipole->Amplitude(0.01, q1, q2);
        //* Cos2Phi_N(r, b2, dipole);
    }
    else    // xH0
    {
        double i0 = gsl_sf_bessel_I0(2.0*b*b2/(l*l));
        double j0 = gsl_sf_bessel_J0(k*r);
        // eq 25
        result = -1.0/(2.0*M_PI) * std::exp(-1.0/(l*l) *(b*b+b2*b2) - r*r/(4.0*l*l));
        
        result *=(((1.0/(l*l) * (b*b + b2*b2) + l*l*k*k - r*r/(4.0*l*l)) * i0 - 2.0*b*b2 /(l*l) * i1 ) * j0  - k*r*i0*j1)
        * dipole->Amplitude(0.01, q1, q2);
        //* Cos2Phi_N(r, b2, dipole);
    }
    //
    
    
    int Nc=3;
    double as=0.3;
    result *= 2.0*Nc/(l*l*l*l*as*M_PI) * b2*r;
    
    return result;
    //return r*result;   // jacobian as we integrate log r
}
