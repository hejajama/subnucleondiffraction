/*
 * Calculate DIS cross section
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 */

#include "dis.hpp"
#include <gsl/gsl_integration.h>
#include "vector.hpp"
#include <amplitudelib/wave_function.hpp>

using Amplitude::SQR;

const int MAXITER_INT = 10;

DIS::DIS(DipoleAmplitude* dipole_)
{
    dipole=dipole_;
}

double DIS::F2(double qsqr, double xbj)
{
    double xs_l, xs_t;
    xs_l = PhotonProtonCrossSection(qsqr, xbj, L);
    xs_t = PhotonProtonCrossSection(qsqr, xbj, T);
    return qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
}

struct Inthelper_dis
{
    double Qsqr, xbj;
    WaveFunction *wf;
    Polarization pol;
    DipoleAmplitude* amp;
    double b;
    double r;
    double theta_b;
};

double inthelperf_dis_b(double b, void* p);
double inthelperf_dis_r(double r, void* p);
double inthelperf_dis_theta_b(double theta_b, void* p);
double inthelperf_dis_theta_r(double theta_r, void* p);

double DIS::PhotonProtonCrossSection(double qsqr, double xbj, Polarization pol)
{
    Inthelper_dis par;
    par.Qsqr=qsqr;
    par.xbj=xbj;
    par.pol = pol;
    par.amp = dipole;
    
    gsl_function fun;
    fun.function = &inthelperf_dis_b;
    fun.params = &par;
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_INT);
    int status = gsl_integration_qag(&fun, 0.000001, 100, 0, 0.01,
                                     MAXITER_INT, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    return result;
}

double inthelperf_dis_b(double b, void* p)
{
    Inthelper_dis* par = (Inthelper_dis*)p;
    par->b = b;
    gsl_function fun;
    fun.params=par;
    fun.function = inthelperf_dis_theta_b;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_INT);
    int status = gsl_integration_qag(&fun, 0.0, 2.0*M_PI, 0, 0.01,
                                     MAXITER_INT, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    return b*result;
}

double inthelperf_dis_theta_b(double theta_b, void* p)
{
    Inthelper_dis* par = (Inthelper_dis*)p;
    par->theta_b = theta_b;
    gsl_function fun;
    fun.params=par;
    fun.function = inthelperf_dis_r;
    
    
    
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_INT);
    int status = gsl_integration_qag(&fun, 0, 999, 0, 0.01,
                                     MAXITER_INT, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    return result;
}



double inthelperf_dis_r(double r, void* p)
{
    Inthelper_dis* par = (Inthelper_dis*)p;
    par->r = r;
    gsl_function fun;
    fun.params=par;
    fun.function = inthelperf_dis_theta_r;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_INT);
    int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, 0, 0.01,
                                     MAXITER_INT, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    return r*result;
}

double inthelperf_dis_theta_r(double theta_r, void* p)
{
    Inthelper_dis* par = (Inthelper_dis*)p;
    Vec b(par->b*std::cos(par->theta_b), par->b*std::sin(par->theta_b));
    Vec r(par->r*std::cos(theta_r), par->r*std::sin(theta_r));
    Vec q1 = r*0.5;
    q1 += b;
    Vec q2 = r*(-0.5);
    q2 += b;
    
    double result = par->amp->Amplitude(par->xbj, q1, q2);
    if (par->pol==L)    // Longitudinal
        result *= par->wf->PsiSqr_L_intz(par->Qsqr, par->r);
    else if (par->pol==T)   // Transverse
        result *= par->wf->PsiSqr_T_intz(par->Qsqr, par->r);
    
    return result;

    
}

DipoleAmplitude* DIS::GetDipole()
{
    return dipole;
}
