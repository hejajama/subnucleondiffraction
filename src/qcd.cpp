/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "qcd.hpp"
#include "subnucleon_config.hpp"
#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>

using namespace std;



QCD::QCD()
{
    // Initialize default values
    nc=3;
    nf=3;
    rc=RUNNING;
    lqcd = 0.241;
    maxalpha=0.7;
}

/*
 * Strong coupling
 */
double QCD::Alphas(double Qsqr)
{
	if (rc == FIXED)
        return 0.2;
	
	if (Qsqr <= lqcd*lqcd)
		return maxalpha;

    double alpha = 12.0*M_PI/( (11.0*Nc()-2.0*Nf())*log(Qsqr/(lqcd*lqcd)) );
    if (alpha > maxalpha)
        return maxalpha;

    return alpha;
}

double QCD::Cf()
{
    return (Nc()*Nc()-1.0)/(2.0*Nc());
}

int QCD::Nc()
{
    return nc;
}

int QCD::Nf()
{
    return nf;
}

void QCD::SetRunningCoupling(RunningAlphas rc_)
{
    rc=rc_;
}


double W_S(double r, int A);
double w_s_normalization=1;
double w_s_normalization_A=-1;

struct inthelper_ta
{
        double b;
        int A;
};

// WS initialization, calculates normalization s.t.
// \int d^3 b WS(b)=1
double inthelperf_ws(double r, void* p)
{
        return r*r*W_S(r, *((int*)p));
}

// Return W_S(\sqrt{ b^2 + z^2 } )
double inthelperf_ta(double z, void* p)
{
        inthelper_ta* par = (inthelper_ta*)p;
        return W_S(std::sqrt( SQR(par->b) + SQR(z) ), par->A);
}



double T_A(double b, int A)
{
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
        double res, abserr;
        inthelper_ta par; par.A=A; par.b=b;
        gsl_function f; f.params=&par; f.function=inthelperf_ta;
        int status = gsl_integration_qag(&f, 0, 99, 0, 0.001, 10, GSL_INTEG_GAUSS61,
                w, &res, &abserr);
        gsl_integration_workspace_free(w);
        if(status)
                cerr << "T_A integration failed at " << LINEINFO <<", result " << res
                        << ", relerr " << std::abs(res-abserr)/res <<", A=" << A << endl;
        return 2.0*res; // 2.0 as we integrate z in [0,\infty]
}


void InitializeWSDistribution(int A)
{
        w_s_normalization = 1;
        w_s_normalization_A = A;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
        double res, abserr;
        gsl_function f; f.params=&A; f.function=inthelperf_ws;
        int status = gsl_integration_qag(&f, 0, 99, 0, 0.001, 10, GSL_INTEG_GAUSS61,
                w, &res, &abserr);
        if (status)
                cerr << "WS normalization integral failed at " << LINEINFO <<", result "
                        << res << " relerr " << std::abs(res-abserr)/res << endl;

        res *= 4.0*M_PI;

        w_s_normalization = 1.0/res;
    gsl_integration_workspace_free(w);
}


double W_S(double r, int A)
{       
        if (A != w_s_normalization_A)
        {       
                cerr << "Woods-Saxon distribution is not initialized for A="<< A << " " << LINEINFO << endl;
                return 0;
        }
        // R = 1.12 fm A^(1/3) - 0.86 fm * A^(-1/3)
        double ra = 1.12 * std::pow(A, 1.0/3.0) - 0.86 * std::pow(A, -1.0/3.0);
        double delta = 0.54 * FMGEV;
        ra *= FMGEV;    // fm => 1/GeV
        
        return w_s_normalization / (std::exp((r - ra)/delta)+1);
}



