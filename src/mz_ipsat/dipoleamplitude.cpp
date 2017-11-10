/*
 * Solve DGLAP equation with given parametrization
 * and evaluate dipole amplitude
 *
 * Uses exactly the same DGLAP solver (LO_evolution_routine.f and
 * alphaS.f) as is used in the fit code
 *
 * H. Mäntysaari and P. Zurita, 2017
 */

#include "dipoleamplitude.hpp"
#include <string>
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace MZ_ipsat;

    // Parameters are: DipoleAmplitude(C, mu0 [GeV], lambda_g, A_g, m_c [GeV]
    // ipsat
    //DipoleAmplitude amplitude(2.146034445992, 1.1, 0.09665075464199, 2.103826220003, 1.351650642298); //
    //amplitude.SetSaturation(true);
 
    // IPnonsat
    //DipoleAmplitude amplitude(4.939286653112, 1.1, -0.009631194037871, 3.058791613883, 1.342035015621);
    //amplitude.SetSaturation(false);'


///////////////////////////////
// DipoleAmplitude class, wrapper to all DGLAP fortran etc


DipoleAmplitude::DipoleAmplitude(double C_, double mu0_, double lambda_g_, double A_g_, double mc_, double mb_, double mt_)
{
    C=C_; mu0=mu0_; lambda_g = lambda_g_; A_g = A_g_; mc = mc_; mb = mb_; mt = mt_;
    saturation = true;
    B_p=4.0;    // All fits are done with fixed B_p=4
    Nc=3;
    mb=4.75;
    mt=175;
    
    // Init alphas(M_Z=91.1876 GeV) = 0.1183
    alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, mc, mb, mt);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    cppdglap = new EvolutionLO(alphas);
}

DipoleAmplitude::~DipoleAmplitude()
{
    delete cppdglap;
}

double DipoleAmplitude::Alphas_xg(double x, double musqr)
{
    double as_xg=0;
    int coupling = 0;
    double As=0;
    double lambdas=0; //singlet
    
    return cppdglap->alphasxG(x, musqr, mu0, 0, A_g, lambda_g, As, lambdas);
}

double DipoleAmplitude::Alphas(double Q)
{
    return alphas->value(Q);
}

double DipoleAmplitude::xg(double x, double musqr)
{
    return cppdglap->xG(x, musqr, mu0, 0, A_g, lambda_g, 0 , 0);
}

double DipoleAmplitude::N(double r, double xbj, double b)
{
    double musqr = mu0*mu0 + C / (r*r);
    double exponent = M_PI*M_PI / (2.0 * Nc) * r*r * Alphas_xg(xbj, musqr) * Tp(b);
    
    if (!saturation)     // IPnonsat
        return exponent;
    else
        return 1.0 - std::exp(-exponent);

}


struct inthelper_bint
{
    double r,x;
    DipoleAmplitude* ipsat;
};
double inthelperf_bint(double b, void* p)
{
    inthelper_bint* par = (inthelper_bint*) p;
    return b*par->ipsat->N(par->r, par->x, b);
}


double DipoleAmplitude::N_bint(double r, double xbj)
{
    double musqr =mu0*mu0 + C / (r*r);
    if (!saturation)
    {
        return 2.0 * M_PI * B_p * N(r, 0, xbj);
    }
    

    
    double a = M_PI*M_PI / (2.0 * Nc) * r*r * Alphas_xg(xbj, musqr)  / (2.0 * M_PI * B_p);
    if (a==0) // Basically so small r that xg =0 as we are outside the dglap evolution grid
        return 0;
    
    gsl_sf_result sinres;
    gsl_sf_result cosres;
    int sinint = gsl_sf_Shi_e(a, &sinres);
    // No overflows
    int cosint = gsl_sf_Chi_e(a, &cosres);
    
    
    if (cosres.val < 1e3 and sinres.val < 1e3 and !cosint and !sinint)
    {
        // Dont trust result if this condition is not true, as we have to compute consint - sinint
        // And these functions grow very rapidly with argument
        // Here we use the analytical result, otherwise we fall back to numerics
        // in the region where the contribution anyway is small
        return 2.0*M_PI*B_p * ( M_EULER - cosres.val + log(a) + sinres.val);
    }
    
    gsl_function fun; fun.function=inthelperf_bint;
    inthelper_bint par;
    par.r=r; par.x=xbj;
    par.ipsat = this;
    fun.params=&par;
    
    double acc = 0.00001;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(500);
    int status = gsl_integration_qag(&fun, 0, 999, 0, acc,
                                     500, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    if (status)
        cerr << "bintegral failed in IPsat::DipoleAmplitude_bit with r=" << r <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return 2.0*M_PI*result; //2pi from angular integral
    
}

double DipoleAmplitude::Tp(double b)
{
    return 1.0/(2.0*M_PI*B_p) * std::exp(-b*b / (2.0*B_p));
}


