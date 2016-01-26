/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "ipsat_proton.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <cmath>
#include <string>
#include <sstream>
#include "subnucleon_config.hpp"

// IPsat 2012
extern "C" {
       double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
     };

int IPSAT12_PAR = 2;    // m_c=1.4 GeV

using std::cout; using std::endl;

const double FMGEV = 5.06778;

void Ipsat_Proton::InitializeTarget()
{
    double smallestdist=999; double largestdist=0;
    double tmpdist;

    
    
    quarks.clear();
    quark_bp.clear();
    
    // Sample 3 quarks
    for (int i=0; i<3; i++)
    {
        // Radius from uniform distribution
        double radius=0;
        double maxr = 30;
        
        if (shape == GAUSSIAN)
        {
            if (B_p < 1e-5)
                radius=0;
            else
            {
                do{
                    radius = gsl_rng_uniform(global_rng) * maxr;
                } while (gsl_rng_uniform(global_rng) > GaussianRadiusDistribution(radius));
            }
            // Sample angle
            double angle = 2.0*M_PI*gsl_rng_uniform(global_rng);
            Vec tmpvec(radius*std::cos(angle), radius*std::sin(angle));
            quarks.push_back(tmpvec);
            quark_bp.push_back(B_q);
        }
        else if (shape == EXPONENTIAL)
        {
            // We have to sample x,y,z separately
            double x,y,z;
            do{
                x = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
                y = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
                z = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
            } while (gsl_rng_uniform(global_rng) > ExponentialDistribution(x,y,z));
            Vec tmpvec(x,y,z);
            quarks.push_back(tmpvec);
            quark_bp.push_back(B_q);
        }
    }
    //exit(1);
    
    
}

Ipsat_Proton::Ipsat_Proton()
{
    saturation=true;
    B_p = 4.0; // GeV^-2
    B_q = B_p/3.0;
    shape = GAUSSIAN;
 
    gdist = new DGLAPDist();
    allocated_gdist = true;
    ipsat = IPSAT12;
    skewedness=false;
}
Ipsat_Proton::Ipsat_Proton(DGLAPDist *gd)
{
    saturation=true;
    B_p = 4.0; // GeV^-2
    B_q = B_p/3.0;
    shape = GAUSSIAN;
    
    gdist = gd;
    allocated_gdist = false;
    ipsat = IPSAT06;
    skewedness=false;

}

Ipsat_Proton::~Ipsat_Proton()
{
    if (allocated_gdist)
        delete gdist;
}




double Ipsat_Proton::Amplitude( double xpom, double q1[2], double q2[2])
{
    // quark transveser coodinates are now nucleons[i].GetX() and GetY()
    Vec q(q1[0], q1[1]);
    Vec qbar(q2[0],q2[1]);
    
    Vec b = q + qbar;
    b = b*0.5;
    
    Vec r = q - qbar;
    
    // Need to calculate sum_i T_p(|b-b_i|), where b is the center of the dipole
    // and b_i is the center of the quark i
    // Note that at this point we do not take into account z, so b in geometric average, not
    // center of mass. Should we change this????
    
    // Currently should correspond to calculation arXiV:1011.1988
    
    double tpsum = 0;
    for (unsigned int i=0; i < quarks.size(); i++)
    {
        Vec projection (quarks[i].GetX(), quarks[i].GetY());
        Vec deltab = b - projection;
        tpsum = tpsum + QuarkThickness(deltab.Len(), i);
    }
    
    if (ipsat == IPSAT06)
        {
            double skew=1.0;
            if (skewedness)
            {
                double skew_lambda = LogDerivative_xg(xpom, r.Len());
                skew = Skewedness(skew_lambda);
            }
            //std::cout << "skew at r=" << r.Len() << " is " << skew << std::endl;
        if (saturation)
            return 1.0 - std::exp( - r.LenSqr() * 1.0/quarks.size() * skew*gdist->Gluedist(xpom, r.LenSqr()) * tpsum  );
        else
            return r.LenSqr() * 1.0/quarks.size() * gdist->Gluedist(xpom, r.LenSqr()) * tpsum  ;
        }
    else if (ipsat == IPSAT12)
    {
        if (!saturation)
        {
            std::cerr << "Nonsat is not defined for ipsat2012" << std::endl;
            exit(1);
        }
        // dipole_amplitude(xBj, r, b, parameterSet) gives amplitude 2(1 - exp(c*T(p)))
        // We have to calculate "gluedist" as in case of ipsat06
        double tmpb=0;  double tmpr = r.Len();
        // par 1: m_c=1.27,   2: m_c=1.4
        double n = dipole_amplitude_(&xpom, &tmpr, &tmpb, &IPSAT12_PAR)/2.0;
        
        double c = std::log(1.0-n);
        double tp = 1.0/(2.0*M_PI*4.0)*std::exp(- tmpb / (2.0*4.0));
        c /= tp;
        
        double skew=1.0; // now c contains xg that is modified by skewedness correction if enabled
        if (skewedness)
        {
            double skew_lambda = LogDerivative_xg(xpom, r.Len());
            if (std::isnan(skew_lambda))
            {
                std::cerr << "skew nan at xpom=" << xpom << " r = " << r.Len() << std::endl;
                skew=1.0;
            }
            else
                skew = Skewedness(skew_lambda);
        }
        
        return 1.0 - std::exp(skew*c * 1.0/quarks.size()*tpsum);
    }
    else
    {
        std::cerr << "UNKNOWN IPSAT VERSION!" << std::endl;
        exit(1);
    }

}


/* extract collinear factorization gluon distribution from ipsat
 */
double Ipsat_Proton::xg(double x, double r)
{
    if (ipsat == IPSAT06)
    {
        double gd = gdist->Gluedist(x, r*r);
        double musqr = 4.0/(r*r)+1.17;
        double as = 12.0*M_PI/( (33.0-2.0*3.0)*log(musqr/(0.2*0.2) ) );
        
        return 2.0*3.0/(M_PI*M_PI*as ) * gd;
    }
    else if (ipsat == IPSAT12)
    {
        double tmpb = 0;
        
        double n = dipole_amplitude_(&x, &r, &tmpb, &IPSAT12_PAR)/2.0;
        double exp = std::log(1.0-n);
        // exp is -pi^2 r^2/(2Nc) as xg T(0)
        double musqr = 4.0/(r*r) + 1.428;   // 1.428 corresponds to m_c=1.4 GeV
        double as = 12.0*M_PI / ( (33.0-2.0*3.0)*log(musqr/(0.156*0.156) ) );
        
        return -2.0*3.0/ (M_PI*M_PI * as* r*r * 1.0/(2.0*M_PI*4.0))*exp;
        
    }
    
}

/*
 * d ln xg / d ln (1/x)
 */
double dhelperf_xg(double y, void* p);
struct dhelper_xg { Ipsat_Proton* proton; double r; };
double Ipsat_Proton::LogDerivative_xg(double x, double r)
{
    gsl_function F;
    F.function=&dhelperf_xg;
    dhelper_xg par; par.proton=this; par.r=r;
    F.params = &par;
    double result,abserr;
    double y = std::log(1.0/x);
    gsl_deriv_central (&F, y, 0.1 , &result, &abserr);
    
    return result;

}

double dhelperf_xg(double y, void* p)
{
    dhelper_xg *par = (dhelper_xg*)p;
    double x = std::exp(-y);
    return std::log(par->proton->xg(x, par->r));
}


std::vector<Vec> &Ipsat_Proton::GetQuarks()
{
    return quarks;
}

std::vector<double> Ipsat_Proton::GetRadii()
{
    std::vector<double> radii;
    for (int i=0; i<quark_bp.size(); i++)
        radii.push_back(std::sqrt(2.0*quark_bp[i]));
    return radii;
}

double Ipsat_Proton::QuarkThickness(double r, int i)
{
    double bp = quark_bp[i];
    return 1.0/(2.0*M_PI*bp)*std::exp(- r*r / (2.0*bp));
}

/*
 * Quark distances from the origin are sampled from this distribution
 */
double Ipsat_Proton::GaussianRadiusDistribution(double r)
{
    return std::sqrt(std::exp(1) / std::exp(B_p))*r*std::exp( - r*r / (2.0*B_p));
}
    // Factor 0.5 should not matter, in just puts the probability always < 1 to accept
    // a sampled radius. Any factor ]0,1] should work (TEST!)

double Ipsat_Proton::ExponentialDistribution(double x, double y, double z)
{
    //a = \sqrt{12}/R_p = 3.87, with R_p = 0.895 from
    //http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1
    double a = B_p;
    
    return std::exp( -a * std::sqrt( x*x + y*y + z*z ) );
    
}



        
std::string Ipsat_Proton::InfoStr()
{
    std::stringstream ss;
    ss << "IPsat proton consists of quarks at coordinates ";
    for (int i=0; i<quarks.size(); i++)
    {
        ss << "(" << quarks[i].GetX() << ", " << quarks[i].GetY() << "), r=" << std::sqrt(2.0*quark_bp[i]) <<", B_q=" << quark_bp[i] ;
    }
    ss << ". Proton radius " << std::sqrt(2.0*B_p) << " GeV^-1, B_p=" << B_p << " ";
    if (shape == GAUSSIAN)
        ss << "Gaussian shape";
    else if (shape == EXPONENTIAL)
        ss << "Exponential shape, a=B_p";
    
    ss << ". Saturation: ";
    if (saturation)
        ss << " enabled";
    else ss << "disabled";
    ss << ". ";
    if (ipsat == IPSAT06)
        ss << "IPsat version: 2006 (KMW)";
    else if (ipsat == IPSAT12)
        ss << "IPsat version: 2012";
    ss << ". Skewedness in dipole amplitude: ";
    if (skewedness)
        ss << " Enabled";
    else
        ss << "Disabled";
    return ss.str();
}


void Ipsat_Proton::SetProtonWidth(double bp)
{
    B_p = bp;
}

void Ipsat_Proton::SetQuarkWidth(double bq)
{
    B_q = bq;
}

void Ipsat_Proton::SetShape(Proton_shape s)
{
    shape = s;
}

double Ipsat_Proton::Amplitude(double xpom, Vec q1, Vec q2)
{
    double quark[2] = {q1.GetX(), q1.GetY() };
    double antiquark[2] = {q2.GetX(), q2.GetY() };
    return Amplitude(xpom, quark, antiquark);
}

double Ipsat_Proton::Skewedness(double lambda)
{
    if (lambda+4.0 > GSL_SF_GAMMA_XMAX or 2.0*lambda+3 > GSL_SF_GAMMA_XMAX)
    {
        std::cerr << "Cant calculate skewedness, overflow, lambda=" << lambda << std::endl;
        return 1.0;
    }
    gsl_sf_result gamma_lambda_52;
    int status1 =gsl_sf_gamma_e(lambda+5.0/2.0, &gamma_lambda_52);
    gsl_sf_result gamma_lambda_4;
    int status2 =gsl_sf_gamma_e(lambda+4.0, &gamma_lambda_4);
    if (status1 or status2)
    {
        std::cerr << "Gamma function evaluation failed, Gamma(" <<lambda+5.0/2.0 << ")=" << gamma_lambda_52.val <<", Gamma(" << lambda+4.0 << ")=" << gamma_lambda_4.val << ", lambda=" << lambda << std::endl;
        return 1.0;
    }
    return std::pow(2.0, 2.0*lambda+3.0)/std::sqrt(M_PI) * gamma_lambda_52.val/gamma_lambda_4.val;
}

void Ipsat_Proton::SetSkewedness(bool s)
{
    skewedness = s;
}