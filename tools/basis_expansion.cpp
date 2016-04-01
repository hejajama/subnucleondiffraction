/*
 * Expand proton density profile in basis functions
 * to see what components are most dominant
 * Basis functions are radial Bessel functions and
 * for the angular part ~ exp(im phi)
 *
 * Ref. for example http://lmb.informatik.uni-freiburg.de/papers/download/wa_report01_08.pdf 
 * Chapter 2
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 */

#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include "../src/dipole.hpp"
#include <string>
#include <sstream>
#include <complex>
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <tools/tools.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_bessel.h>
using namespace std;
using Amplitude::SQR;
gsl_rng* global_rng;

const double MAXR = 20;
const int INTERVALS = 20;
const double INTACCURACY = 0.001;
const double INTACCURACY_ABS = 1e-6;
const int CONFIGURATIONS = 1000;

const unsigned int MAXN = 10;
const unsigned int MAXM = 10;

// Calculate integral (41) form Ref.
std::complex<double> ExpansionCoefficient(unsigned int n,  unsigned int m, DipoleAmplitude* amp);
std::complex<double> BasisFunction(unsigned int n, unsigned int m, double r, double phi);

int main(int argc, char* argv[])
{
    //gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    DipoleAmplitude *amp;

    // Arguments similarly as in main.cpp for -dipole
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-dipole")
        {
            if (string(argv[i+1])=="ipsatproton")
            {
                amp = new Ipsat_Proton;
                ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+2]));
                ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+3]));
                ((Ipsat_Proton*)amp)->SetShape(GAUSSIAN);
                if (argc > i+4)
                {
                    if (string(argv[i+4])=="fluxtube")
                    {
                        ((Ipsat_Proton*)amp)->SetStructure(CENTER_TUBES);
                        ((Ipsat_Proton*)amp)->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                    }
                    else if (string(argv[i+4]).substr(0,1)!="-")
                    {
                        cerr << "Unknown ipsatproton option " << argv[i+4] << endl;
                        exit(1);
                    }
                }
            }
            else if (string(argv[i+1])=="ipglasma")
                amp = new IPGlasma(argv[i+2]);
        }
    }
    
    //amp->InitializeTarget();
    
    // Calculate coefficients
    cout << "n  m  P_nm" << endl;
    for (int n=1; n<MAXN; n+=1)
    {
        for (int m=0; m<=MAXM; m+=1)
        {
            double values_real[CONFIGURATIONS];
            double values_imag[CONFIGURATIONS];
	    double values_abs[CONFIGURATIONS];
            for (int conf=0; conf < CONFIGURATIONS; conf++)
            {
                gsl_rng_set(global_rng, conf);
                amp->InitializeTarget();
                complex<double> pnm = ExpansionCoefficient(n, m, amp);
                values_real[conf]=pnm.real();
                values_imag[conf]=pnm.imag();
		values_abs[conf]=norm(pnm); 
            }
            double real = gsl_stats_mean( values_real, 1, CONFIGURATIONS );
            double imag = gsl_stats_mean( values_imag, 1, CONFIGURATIONS );
            double abs = gsl_stats_mean( values_abs, 1, CONFIGURATIONS );
            cout << n << " " << m << " " << abs << " " << real << " " << imag << endl;
        }
    }
    
    delete amp;
    

}

struct inthelper_coef { unsigned int n; unsigned int m; double r; DipoleAmplitude* amp; bool imag; };
double inthelperf_r(double r, void* p);
double integrand(double phi, void* p);

std::complex<double> ExpansionCoefficient(unsigned int n, unsigned int m, DipoleAmplitude* amp)
{
    inthelper_coef par;
    par.n=n;
    par.m=m;
    par.amp=amp;
    
    
    double real,imag;
    
    gsl_function fun;
    fun.function = &inthelperf_r;
    fun.params = &par;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    
    // Real part
    double err_real;
    par.imag = false;
    int status = gsl_integration_qag(&fun, 0, MAXR, INTACCURACY_ABS, INTACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &real, &err_real);
    // Imag
    double err_imag;
    par.imag = true;
     status = gsl_integration_qag(&fun, 0, MAXR, INTACCURACY_ABS, INTACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &imag, &err_imag);
    
    gsl_integration_workspace_free(w);
    
    
    complex<double> pnm(real,imag);
    
    return pnm;
    
}

double inthelperf_r(double r, void* p)
{
    inthelper_coef *par = (inthelper_coef*)p;
    par->r = r;
    
    gsl_function fun;
    fun.function = &integrand;
    fun.params = par;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    
    double result, err;
    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, INTACCURACY_ABS, INTACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &err);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double integrand(double phi, void* p)
{
    inthelper_coef *par = (inthelper_coef*)p;
    
    complex<double> basis_c = BasisFunction(par->n, par->m, par->r, phi);
    double basis;
    if (par->imag)
        basis = basis_c.imag();
    else
        basis = basis_c.real();
    
    Vec b(par->r * cos(phi), par->r * sin(phi));
    
    return par->r * basis * par->amp->Density(b);
}

std::complex<double> BasisFunction(unsigned int n, unsigned int m, double r, double phi)
{
    // Angular part
    std::complex<double> exponent(0, m*phi);
    std::complex<double> angular = 1.0/sqrt(2.0*M_PI) * exp(exponent);
    double xmn = gsl_sf_bessel_zero_Jnu(m, n);
    
    
    // Radial part bessel function
    double k= xmn / MAXR;
    gsl_sf_result res;
    int status = gsl_sf_bessel_J0_e(k*r, &res);
    if (status)
    {
        cerr << "GSL bessel error at r=" << r << endl;
        exit(1);
    }
    double bessel = res.val;
    
    // Use zero value boundary condition
    double normalization_nm = SQR(MAXR)/2.0 * SQR(gsl_sf_bessel_Jn(m+1,xmn));
    double radial = 1.0 / sqrt( normalization_nm ) * gsl_sf_bessel_Jn(m, k*r);
    
    return radial*angular;
}
