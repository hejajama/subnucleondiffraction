/*
 * Overlap between the photon and the vector meson wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "gaus_lc.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

const REAL ZINTACCURACY=0.001;
const int MAXITER_ZINT=1000;
const REAL MINZ=0.00001;  // Integration limits
const REAL MAXZ=0.9999;

using Amplitude::SQR;
const int NC=3;

using Amplitude::ALPHA_e;

//const double ALPHA_e = 1.0/137.035999679;
const double e = sqrt(4.0*M_PI*ALPHA_e);

GausLC::GausLC(REAL e_f_, REAL N_T_, REAL N_L_, REAL R_T_, REAL R_L_, 
    REAL m_f_, REAL M_V_, int delta_=0 )
{
    e_f=e_f_;
    m_f=m_f_; M_V=M_V_;
    N_T=N_T_; N_L=N_L_;
    R_T=R_T_; R_L=R_L_;
    delta=delta_; 
     
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
GausLC::GausLC(std::string file)
{
    std::ifstream f(file.c_str());
  
    if (!f.is_open())
    {
        std::cerr << "Could not open file " << file << "!";
        return;
    }
    std::string line;
    
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
        if (line.substr(0,3)=="R_T") 
            R_T=StrToReal(line.substr(4));
        if (line.substr(0,3)=="R_L") 
            R_L=StrToReal(line.substr(4));
        if (line.substr(0,3)=="del") 
            delta=StrToInt(line.substr(4));
    }
    f.close();

}


/*
 * Transversially polarized component of the overlap
 */

REAL GausLC::PsiSqr_T(REAL Qsqr, REAL r, REAL z)
{
    REAL result;
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    result = SQR(m_f)*gsl_sf_bessel_K0(epstmp*r)*Psi_T(r,z) 
        - (SQR(z)+SQR(1-z))*epstmp*gsl_sf_bessel_K1(epstmp*r)*Psi_T_DR(r,z);
    result *= e_f*e*NC/(M_PI*z*(1-z));
    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

REAL GausLC::PsiSqr_L(REAL Qsqr, REAL r, REAL z)
{
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    REAL result;
    result = M_V*Psi_L(r,z)+delta/(M_V*z*(1-z))*(SQR(m_f)*Psi_L(r,z) 
        - 1/r*Psi_L_DR(r,z) - Psi_L_D2R(r,z));
    result *= e_f*e*NC/M_PI*2*sqrt(Qsqr)*z*(1-z)*gsl_sf_bessel_K0(epstmp*r);
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
struct zinthelper{
    GausLC *vm_p;
    double  Qsqr;
    double  r;
};

REAL  zhelperfunc_lc_T(REAL z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_T(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

REAL  zhelperfunc_lc_L(REAL z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_L(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

REAL GausLC::PsiSqr_T_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfunc_lc_T;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, MINZ, MAXZ,  0, ZINTACCURACY, 
        &result, &abserr, &eval);
    //gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    //status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
    //    MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    //gsl_integration_workspace_free(ws);
    
    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << status << " (transverse, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ")" << std::endl;}

    result*=1.0/(4.0*M_PI); // Normalization
    return result;
}

REAL GausLC::PsiSqr_L_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfunc_lc_L;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY, 
        &result, &abserr, &eval);
    //gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    //int status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
    //    MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    //gsl_integration_workspace_free(ws);
    
    if(status){ std::cerr<< "\\int z in GausLC failed: code " 
        << status << " (longitudinal, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ")" << std::endl;}

    result*=1.0/(4.0*M_PI); // Normalization
    return result;
}

/*
 * Gaussian lightcone functions 
 */
REAL GausLC::Psi_T(REAL r, REAL z)
{
    return N_T*SQR(z*(1.0-z))*exp(-SQR(r)/(2.0*SQR(R_T)));
}

REAL GausLC::Psi_L(REAL r, REAL z)
{
    return N_L*z*(1.0-z)*exp(-SQR(r)/(2.0*SQR(R_L)));
}

// \partial_r Psi_T(r,z)
REAL GausLC::Psi_T_DR(REAL r, REAL z)
{
    return -1.0/(SQR(R_T))*N_T*r*SQR(z*(1-z))*exp(-SQR(r)/(2.0*SQR(R_T)));
}

// \partial_r PSI_L(r,z)
REAL GausLC::Psi_L_DR(REAL r, REAL z)
{
    return -1.0*r/(SQR(R_L))*N_L*z*(1-z)*exp(-SQR(r)/(2.0*SQR(R_L)));
}

// \partial^2_r PSI_L(r,z)
REAL GausLC::Psi_L_D2R(REAL r, REAL z)
{
    return N_L*(z-1)*z/pow(R_L,4) * (SQR(r)-SQR(R_L)) * exp(-SQR(r)/(2.0*SQR(R_L)));
}

REAL GausLC::MesonMass()
{
    return M_V;
}

std::string GausLC::GetParamString()
{
    std::stringstream str;
    str << "e_f = "
    << e_f << ", m_f = " << m_f << ", M_V = " << M_V << ", N_T = "
    << N_T << ", N_L = " << N_L << ", R_T = " << R_T << ", R_L = " << R_L
    << ", delta = " << delta << " ";
    return str.str();
}

std::ostream& operator<<(std::ostream& os, GausLC& ic)
{
    return os << " Wave function overlap, Gaus-LC, Params: "
        << ic.GetParamString() << " .";
        
}

