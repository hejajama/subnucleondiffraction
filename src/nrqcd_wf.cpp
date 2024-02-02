/*
 * Overlap between the photon and the vector meson wave functions
 * NRQCD parametrization
 * T. Lappi, H. Mäntysaari, J. Penttala, arXiv:2006.02830
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2021
 */
 
#include "nrqcd_wf.hpp"
#include "subnucleon_config.hpp"
#include "qcd.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "subnucleon_config.hpp"


inline double ABS(double x) { if (x<0) return -x; else return x; }
const int NC=3;

const double ZINTACCURACY=ZINT_RELACCURACY;
const int MAXITER_ZINT=ZINT_INTERVALS;

const double e = sqrt(4.0*M_PI*ALPHA_e);


NRQCD_WF::NRQCD_WF(double A_, double B_, double ef_, double mc_, double MV_)
{
    mc=mc_;MV=MV_;A=A_;B=B_; ef=ef_;
    
}


// z dep overlaps not defined, as ~delta(z-1/2)
double NRQCD_WF::PsiSqr_T(double Qsqr, double r, double z)
{
    std::cerr << "NRQCD WF does not support z dependent overlap!" << std::endl;
    return 0;
}

/*
 * Longitudinally polarized component of the overlap
 */

double NRQCD_WF::PsiSqr_L(double Qsqr, double r, double z)
{
    std::cerr << "NRQCD WF does not support z dependent overlap!" << std::endl;
    return 0;
}


/*
 * z integrated overlaps
 * \int dz/(4pi) Psi^*Psi * exp(i(1/2-z)*r*Delta)
 * Note: does not depend on the overall sign in the exponent
 */

double NRQCD_WF::PsiSqr_T_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
{
    double result=0;
    double eps = std::sqrt(mc*mc+Qsqr/4.);
    
    double K0 =gsl_sf_bessel_K0(eps*r);
    
    result = A*e*ef*std::sqrt(mc*NC)*K0/(2.0*std::sqrt(2.0)*M_PI);
    
    // Rel correction
    result += B / (16.0*std::sqrt(2)*std::pow(mc,3./2.)*M_PI*2.*eps) * e*ef*
    (
        -8*std::sqrt(NC)*(2.0*mc*mc+Qsqr)*r*gsl_sf_bessel_K1(eps*r)
     + std::sqrt(NC)*2.*eps*K0*(28.0+r*r*(8.0*mc*mc+Delta*Delta) + r*r*Delta*Delta*std::cos(2.0*phi_r_Delta))
     );
    
    
    return -result; // Overall sign included here, but it does not matter in VM production
}

double NRQCD_WF::PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
{
    double result=0;
    double eps = std::sqrt(mc*mc+Qsqr/4.);
    
    double K0 =gsl_sf_bessel_K0(eps*r);
    
    result = A*e*ef*std::sqrt(NC*Qsqr)*K0/(8.0*std::sqrt(2.0*mc)*M_PI);
    
    result += B/(64*std::sqrt(2)*std::pow(mc,5./2.)*M_PI) * e*ef*std::sqrt(NC*Qsqr)*(
            -4.0*Qsqr*r*gsl_sf_bessel_K1(eps*r)/std::sqrt(4.0*mc*mc+Qsqr)
            + K0*(36.0+r*r*(8.0*mc*mc+Delta*Delta) +r*r*Delta*Delta*std::cos(2.0*phi_r_Delta))
    );
    
    return -2.0*result; // 2.0 from helicity sum
}



double NRQCD_WF::MesonMass()
{
    return MV;
}

std::string NRQCD_WF::GetParamString()
{
    std::stringstream str;
    str << "A = "
    << A << " GeV^(3/2), B = " << B << " GeV^(7/2), M_V = " << MV << " GeV, m_c = "
    << mc << " GeV";
    return str.str();
}

std::ostream& operator<<(std::ostream& os, NRQCD_WF& ic)
{
    return os << "Wave function overlap: NRQCD, Params: "
        << ic.GetParamString() << " .";
        
}

