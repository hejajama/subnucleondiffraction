/*
 * Virtual photon wave function
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2014
 */
 
#include "dvcs_photon.hpp"
#include "subnucleon_config.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#define LINEINFO __FILE__ << ":" << __LINE__

const double ZINTACCURACY=0.0001;
const int MAXITER_ZINT=700;
const double MINZ=0.000001;  // Integration limits
const double MAXZ=0.99999;

using namespace std;


DVCSPhoton::DVCSPhoton()
{
	// Initialize with light quarks
	
    // Quark charges, 0=u, 1=d, 2=s
    // Parameters same as in Ref. 0902.1112
    e_f.push_back(2.0/3.0); e_f.push_back(-1.0/3.0); e_f.push_back(-1.0/3.0);
    m_f.push_back(0.14); m_f.push_back(0.14); m_f.push_back(0.14);
    
}

/*
 * Transversially polarized component of the overlap
 */

double DVCSPhoton::PsiSqr_T(double Qsqr, double r, double z)
{
    double result=0;
    for (unsigned int f=0; f<e_f.size(); f++)     // Sum quar flavors
    {
        double epstmp=Epsilon(Qsqr,z,f);
        result += SQR(e_f[f])*(
            (SQR(z)+SQR(1.0-z))*epstmp
                    *gsl_sf_bessel_K1(epstmp*r)*m_f[f]*gsl_sf_bessel_K1(m_f[f]*r)
                + SQR(m_f[f]*gsl_sf_bessel_K0(epstmp*r))
            );
    }
    result *= 2.0*qcd.Nc()/M_PI*ALPHA_e;
    
    if (isnan(result))
    {
        std::cerr << "PsiSqr_T=NaN at r=" << r <<", Q^2=" << Qsqr <<", z=" << z << " with quarks: " << endl;
        for (unsigned int f=0; f<e_f.size(); f++)
        {
            std::cerr << "m=" << m_f[f] << ", e_f=" << e_f[f] << std::endl;
            exit(1);
        }
    }

    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

double DVCSPhoton::PsiSqr_L(double Qsqr, double r, double z)
{
    return 0;


}

/* 
 * Wave function overlap integrated over z=[0,1]
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct DVCSzinthelper{
    DVCSPhoton *vm_p;
    double  Qsqr;
    double  r;
};

double  DVCSzhelperfuncT(double z, void * p){
  return ((DVCSzinthelper*)p)->vm_p->PsiSqr_T(
							((DVCSzinthelper*)p)->Qsqr,
							((DVCSzinthelper*)p)->r,
							z);
}

double DVCSPhoton::PsiSqr_L_intz(double Qsqr, double r)
{
    return 0;
}

double DVCSPhoton::PsiSqr_T_intz(double Qsqr, double r)
{
    double result,abserr;
    DVCSzinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&DVCSzhelperfuncT;
    int_helper.params=&zintpar;
    
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, MINZ, MAXZ , 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    if (isnan(result))
    {
        std::cerr << "z integral in Photon = NaN, Qsqr=" << Qsqr << ", r=" << r << " at " << LINEINFO << std::endl;
        exit(1);
    }

    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << status << " result " << result << " (transverse, Qsqr=" << Qsqr << ", r=" << r
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}
  
    return result;
}




double DVCSPhoton::Epsilon(double Qsqr, double z, int f)
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

void DVCSPhoton::SetQuark(Parton p, double mass)
{
	// Clear 
	e_f.clear();
	m_f.clear();
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
    if (m_f[0] < 0 or m_f[0] > 1000 or e_f[0] < -1 or e_f[0]>1)
    {
        std::cerr << "Quark m=" << m_f[0] << ", frac.charge=" << e_f[0] << " was provided to VirtualPhoton, this is crazy! " << LINEINFO << std::endl;
        exit(1);
    }

}

std::string DVCSPhoton::GetParamString()
{
    std::stringstream str;
    for (unsigned int f=0; f<e_f.size(); f++)
		str << "e_f[" << f << "]=" << e_f[f]  << ", m_f[" << f <<"]=" << m_f[f] << " ";
    return str.str();
}

std::ostream& operator<<(std::ostream& os, DVCSPhoton& ic)
{
    return os << " DVCS wave function. Params: "
        << ic.GetParamString() << " .";
        
}

