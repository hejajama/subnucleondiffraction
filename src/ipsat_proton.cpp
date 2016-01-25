/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "ipsat_proton.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <string>
#include <sstream>
#include "subnucleon_config.hpp"

// IPsat 2012
extern "C" {
       double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
     };

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
    ipsat = IPSAT06;
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
        if (saturation)
            return 1.0 - std::exp( - r.LenSqr() * 1.0/quarks.size() * gdist->Gluedist(xpom, r.LenSqr()) * tpsum  );
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
        double tmpb=0; int par=2; double tmpr = r.Len();
        // par 1: m_c=1.27,   2: m_c=1.4
        double n = dipole_amplitude_(&xpom, &tmpr, &tmpb, &par)/2.0;    // last number 1: m_c=1.27, 2: 1.4
        double c = std::log(1.0-n);
        double tp = 1.0/(2.0*M_PI*4.0)*std::exp(- tmpb / (2.0*4.0));
        c /= tp;
        return 1.0 - std::exp(c * 1.0/quarks.size()*tpsum);
    }
    else
    {
        std::cerr << "UNKNOWN IPSAT VERSION!" << endl;
        exit(1);
    }

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