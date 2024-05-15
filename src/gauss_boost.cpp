/*
 * Overlap between the photon and the vector meson wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010, 2015
 */
 
#include "gauss_boost.hpp"
#include "subnucleon_config.hpp"
#include <gsl/gsl_math.h>
#include "qcd.hpp"
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

const REAL ZINTACCURACY=ZINT_RELACCURACY;
const int MAXITER_ZINT=ZINT_INTERVALS;


//const double ALPHA_e = 1.0/137.035999679;
const double e = sqrt(4.0*M_PI*ALPHA_e);


BoostedGauss::BoostedGauss(REAL e_f_, REAL N_T_, REAL N_L_, REAL R_, 
    REAL m_f_, REAL M_V_, int delta_=1 )
{
    e_f=e_f_;
    m_f=m_f_; M_V=M_V_;
    N_T=N_T_; N_L=N_L_;
    R=R_;
    delta=delta_; 
    MINZ=0.00001;  // Integration limits
    MAXZ=0.9999;
}

/* Constructor to read parameters from a file
 * The file should be written using the following syntax, this function
 * doesn't check that the syntax is correct!
 *
 * Lines starting with # are comments, and parameters are written in units of
 * GeV^n using the syntax
 * key:value
 * Keys are: e_f, m_f, M_V, N_T, N_L, R_T, R_L and del
 */
BoostedGauss::BoostedGauss(std::string file)
{
    std::ifstream f(file.c_str());
  
    if (!f.is_open())
    {
        std::cerr << "Could not open file " << file << " " << LINEINFO << "!" << std::endl ;
        exit(1);
        return;
    }
    std::string line;
    alpha2=0;

    while(!f.eof() )
    {
        std::getline(f, line);
        if (line[0]=='#')   // Comment line
            continue;
        if (line.substr(0,3)=="e_f")
            e_f=StrToReal(line.substr(4));
        if (line.substr(0,3)=="m_f") 
            m_f=StrToReal(line.substr(4));
        if (line.substr(0,3)=="M_V") 
            M_V=StrToReal(line.substr(4));
        if (line.substr(0,3)=="N_T") 
            N_T=StrToReal(line.substr(4));
        if (line.substr(0,3)=="N_L") 
            N_L=StrToReal(line.substr(4));
        if (line.substr(0,1)=="R") 
            R=StrToReal(line.substr(2));
        if (line.substr(0,3)=="del") 
            delta=StrToInt(line.substr(4));
        if (line.substr(0,2)=="S:")
            S=StrToInt(line.substr(2));
        if (line.substr(0,6)=="alpha2")
            alpha2=StrToReal(line.substr(7));      
        else if (line.substr(0,5)=="alpha") // else if, as this would be true for alpha2 as well
            alpha=StrToReal(line.substr(6));
              
    }
    f.close();

    MINZ=0.00001;  // Integration limits
    MAXZ=0.9999;
}


/*
 * Transversially polarized component of the overlap
 */

REAL BoostedGauss::PsiSqr_T(REAL Qsqr, REAL r, REAL z)
{
    REAL result;
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    result = SQR(m_f)*gsl_sf_bessel_K0(epstmp*r)*Psi_T(r,z) 
        - (SQR(z)+SQR(1.0-z))*epstmp*gsl_sf_bessel_K1(epstmp*r)*Psi_T_DR(r,z);
    result *= e_f*e*NC/(M_PI*z*(1.0-z));
    
    //result *= exp(-r*r / (5.068*5.068*1.5*1.5));
    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

REAL BoostedGauss::PsiSqr_L(REAL Qsqr, REAL r, REAL z)
{
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    REAL result;
    result = M_V*Psi_L(r,z)+delta/(M_V*z*(1.0-z))*(SQR(m_f)*Psi_L(r,z)
        - 1.0/r*Psi_L_DR(r,z) - Psi_L_D2R(r,z));
    result *= e_f*e*NC/M_PI*2*sqrt(Qsqr)*z*(1.0-z)*gsl_sf_bessel_K0(epstmp*r);
//result *= exp(-r*r / (5.068*5.068*1.5*1.5));
    return result;
}

/* 
 * Wave function overlap integrated over z=[0,1]
 * Normalization factor 1/(4*Pi) is included here!
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct boostzinthelper{
    BoostedGauss *vm_p;
    double  Qsqr;
    double  r;
};

REAL  boostzhelperfuncT(REAL z, void * p){
  return ((boostzinthelper*)p)->vm_p->PsiSqr_T(
							((boostzinthelper*)p)->Qsqr,
							((boostzinthelper*)p)->r,
							z);
}

REAL  boostzhelperfuncL(REAL z, void * p){
  return ((boostzinthelper*)p)->vm_p->PsiSqr_L(
							((boostzinthelper*)p)->Qsqr,
							((boostzinthelper*)p)->r,
							z);
}

REAL BoostedGauss::PsiSqr_T_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    boostzinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&boostzhelperfuncT;
    int_helper.params=&zintpar;
    
    //int status = gsl_integration_qng(&int_helper, MINZ, MAXZ,  0, ZINTACCURACY, 
    //    &result, &abserr, &eval);
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    if (status and ABS(result)>0.000001) std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr <<
            std::endl;

    result*=1.0/(4.0*M_PI); // Normalization
    return result;
}

REAL BoostedGauss::PsiSqr_L_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    boostzinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&boostzhelperfuncL;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, MINZ, MAXZ, ZINTACCURACY, ZINTACCURACY, 
        &result, &abserr, &eval);
    //gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    //int status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
    //    MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    //gsl_integration_workspace_free(ws);
    
    if (status and ABS(result)>0.0000001) std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr <<
            std::endl;

    result*=1.0/(4.0*M_PI); // Normalization
//    result *= exp(-r*r / (5.068*5.068*1.5*1.5));
    return result;
}

/*
 * Boosted Gaussian functions 
 */
REAL BoostedGauss::Psi_T(REAL r, REAL z)
{
    return N_T*Psi(r,z);
}             

REAL BoostedGauss::Psi_L(REAL r, REAL z)
{
    return N_L*Psi(r,z);
}

// Scalar part of the wave function without normalizatin factor
REAL BoostedGauss::Psi(REAL r, REAL z)
{
    double res = z*(1.0-z)*exp(-SQR(m_f)*SQR(R)/(8.0*z*(1.0-z))
        - 2.0*z*(1-z)*SQR(r)/SQR(R) + SQR(m_f)*SQR(R)/2.0 );
    if (S==1)
        return res;
    
    double g = 2.0 - SQR(m_f*R) 
        + SQR(m_f*R)/(4.*z*(1.-z)) - 4.*z*(1.-z)*SQR(r)/SQR(R);

    if (S ==2)
        return res *(1 + alpha * g);
    else if (S==3)
        return res*(1.+alpha*g+alpha2*(SQR(g) + 4.*(1.-4.*z*(1.-z)*SQR(r)/SQR(R))));
    else
    {
        std::cerr << "Unknown state " << S << " in BoostedGauss::Psi_T" << std::endl;
        return 0;
    }
}

// \partial_r Psi(r,z)
REAL BoostedGauss::Psi_DR(REAL r, REAL z)
{
    if (S==1)
        return -4.0*z*(1.0-z)*r/SQR(R)*Psi(r,z); ;
    
    // We need the "S=1" wave function denoted by G_nS in 0905.0102
    double Gn = z*(1.0-z)*exp(-SQR(m_f)*SQR(R)/(8.0*z*(1.0-z))
        - 2.0*z*(1-z)*SQR(r)/SQR(R) + SQR(m_f)*SQR(R)/2.0 ); 
    double Gn_der = -4.0*z*(1.0-z)*r/SQR(R)*Gn;

    if (S>1)
    {
        double der1 = -8*r*(1.-z)*z*alpha*Gn/SQR(R) 
            + (2.-SQR(m_f*R)+SQR(m_f*R)/(4.*z*(1.-z))-4.*SQR(r)*z*(1.-z)/SQR(R))*alpha*Gn_der;
        if (S==2)
            return Gn_der+der1;
 
        double der2 = (-((32*r*(1-z)*z)/SQR(R)) 
            - (16*r*(1-z)*z*(2-SQR(m_f)*SQR(R) + (SQR(m_f)*SQR(R))/(4*(1-z)*z) - (4*r*r*(1-z)*z)/SQR(R)))/SQR(R)) 
                * alpha2 * Gn;
        double der2_2 = (4 * (1 - (4 * SQR(r) * (1 - z) * z) / SQR(R)) + (2 - SQR(m_f) * SQR(R) + (SQR(m_f) * SQR(R)) / (4 * (1 - z) * z) - (4 * r * r * (1 - z) * z) / SQR(R)) * (2 - SQR(m_f) * SQR(R) + (SQR(m_f) * SQR(R)) / (4 * (1 - z) * z) - (4 * r * r * (1 - z) * z) / SQR(R))) * alpha2 * Gn_der;
        return Gn_der + der1 + der2 + der2_2;
    }
    std::cerr << "Unknown state " << S << " in BoostedGauss::Psi_DR" << std::endl;
    return 0;
}

// \partial_r Psi_T(r,z)
REAL BoostedGauss::Psi_T_DR(REAL r, REAL z)
{
    return N_T*Psi_DR(r,z);
}

// \partial_r PSI_L(r,z)
REAL BoostedGauss::Psi_L_DR(REAL r, REAL z)
{
    return N_L*Psi_DR(r,z);
}

// \partial^2_r PSI_L(r,z)
REAL BoostedGauss::Psi_L_D2R(REAL r, REAL z)
{
    if (S>2)
    {
        std::cerr << "BoostedGauss::Psi_L_D2R not implemented for S>2" << std::endl;
        exit(1);
    }
    return -4.0*z*(1.0-z)/SQR(R)*Psi_L(r,z) 
        + SQR(4.0*z*(1.0-z)*r/SQR(R))*Psi_L(r,z)
            - 8.0*alpha*N_T*SQR(z*(1.0-z))/std::pow(R,4.0)
                * (8.0*SQR(r)*z*(-1.0+z) + SQR(R) )
                * std::exp( 2.0*z*(-1.0+z)*SQR(r/R) + SQR(1.0-2.0*z)*SQR(m_f*R)/(8.0*z*(-1.0+z)) );
    
}

REAL BoostedGauss::MesonMass()
{
    return M_V;
}

std::string BoostedGauss::GetParamString()
{
    std::stringstream str;
    str << "e_f = "
    << e_f << ", m_f = " << m_f << ", M_V = " << M_V << ", N_T = "
    << N_T << ", N_L = " << N_L << ", R = " << R 
    << ", delta = " << delta << " " << " S = " << S << ", alpha = " << alpha
    << ", alpha_2=" << alpha2 ;
    return str.str();
}

std::ostream& operator<<(std::ostream& os, BoostedGauss& ic)
{
    return os << "Wave function overlap: Boosted Gaussian, Params: "
        << ic.GetParamString() << " .";
        
}

