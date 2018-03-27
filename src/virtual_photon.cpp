/*
 * Virtual photon wave function
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2016
 */
 
#include "virtual_photon.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace Amplitude;

const int Nc=3;

const double ZINTACCURACY=0.0001;
const int MAXITER_ZINT=500;
const double MINZ=1e-4;  // Integration limits
const double MAXZ=1.0-MINZ;

VirtualPhoton::VirtualPhoton()
{
	// Initialize with light quarks
	// Parameters from IPsat2012 fit
    // Quark charges, 0=u, 1=d, 2=s, 3=c
    e_f.push_back(2.0/3.0); e_f.push_back(-1.0/3.0); e_f.push_back(-1.0/3.0); e_f.push_back(2.0/3.0);
    m_f.push_back(0.01); m_f.push_back(0.01); m_f.push_back(0.01); m_f.push_back(1.4);
    
}

/*
 * Transversially polarized component of the overlap
 */

//const double RSQR_QSQR_LIMIT = 200000;

double VirtualPhoton::PsiSqr_T(double Qsqr, double r, double z)
{
    //if (r*r > RSQR_QSQR_LIMIT / Qsqr) return 0;
    double result=0;
    for (unsigned int f=0; f<e_f.size(); f++)     // Sum quar flavors
    {
        double epstmp=Epsilon(Qsqr,z,f);
        result += SQR(e_f[f])*(
            (SQR(z)+SQR(1.0-z))*SQR(epstmp)
                    *SQR(gsl_sf_bessel_K1(epstmp*r))
                + SQR(m_f[f]*gsl_sf_bessel_K0(epstmp*r))
            );
    }
    result *= Nc/(2.0*SQR(M_PI))*ALPHA_e;

    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

double VirtualPhoton::PsiSqr_L(double Qsqr, double r, double z)
{
    //if (r*r > RSQR_QSQR_LIMIT / Qsqr) return 0;
    double result=0;
    for (unsigned int f=0; f<e_f.size(); f++)     // Sum quar flavors
    {
        double epstmp=Epsilon(Qsqr,z,f);
        result += SQR(e_f[f])* SQR( gsl_sf_bessel_K0(epstmp*r) );
    }
    result *= 2.0*Nc/SQR(M_PI)*ALPHA_e*Qsqr*SQR(z)*SQR(1.0-z);

    return result;
}

/* 
 * Wave function overlap integrated over z=[0,1]
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct zinthelper{
    VirtualPhoton *vm_p;
    double  Qsqr;
    double  r;
};

double  zhelperfuncT(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_T(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

double  zhelperfuncL(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_L(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

double VirtualPhoton::PsiSqr_T_intz(double Qsqr, double r)
{
    double result,abserr;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncT;
    int_helper.params=&zintpar;
    
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);

    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << status << " (transverse, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}
  
    return result;
}

double VirtualPhoton::PsiSqr_L_intz(double Qsqr, double r)
{
    double result,abserr;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncL;
    int_helper.params=&zintpar;
    
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    if(status){ std::cerr<< "z integral in VirtualPhoton failed: code " 
        << status << " (longitudinal, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}

    return result;
}



double VirtualPhoton::Epsilon(double Qsqr, double z, int f)
{
    return std::sqrt( z*(1.0-z)*Qsqr+SQR(m_f[f]) );
}


/*
 * Change wave function to the specific quark, supported quarks are u,d,s,c,b
 * Note: Only this quark is used once this function is called
 * 
 * A good way to check that everything works is to compute e.g. F2 using the
 * default quark content (u,d,s) and these quarks separately
 * 
 * If mass<0 (default value), then use standard literature value for the masses
 */

void VirtualPhoton::SetQuark(Parton p, double mass)
{
	// Clear 
	e_f.clear();
	m_f.clear();
	
    AddQuark(p, mass);
}

void VirtualPhoton::AddQuark(Parton p, double mass)
{
    double m;
    switch(p)
    {
        case LIGHT:
            m=0.14;
            if (mass>=0)
                m=mass;
            m_f.push_back(m);
            e_f.push_back(2.0/3.0);
            m_f.push_back(m);
            e_f.push_back(-1.0/3.0);
            m_f.push_back(m);
            e_f.push_back(-1.0/3.0);
            break;
        case U:
            m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
            e_f.push_back(2.0/3.0);
            break;
        case D:
            m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
            e_f.push_back(-1.0/3.0);
            break;
        case S:
            m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
            e_f.push_back(-1.0/3.0);
            break;
        case C:
            m_f.push_back(1.27);
            if (mass>=0) m_f[0]=mass;
            e_f.push_back(2.0/3.0);
            break;
        case B:
            m_f.push_back(4.2);
            if (mass>=0) m_f[0]=mass;
            e_f.push_back(-1.0/3.0);
            break;
        default:
            cerr << "Unknown parton " << p << " at " << LINEINFO << endl;
    }
}

std::string VirtualPhoton::GetParamString()
{
    std::stringstream str;
    for (unsigned int f=0; f<e_f.size(); f++)
		str << "e_f[" << f << "]=" << e_f[f]  << ", m_f[" << f <<"]=" << m_f[f] << " ";
    return str.str();
}

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic)
{
    return os << " Virtual photon wave function. Params: "
        << ic.GetParamString() << " .";
        
}

