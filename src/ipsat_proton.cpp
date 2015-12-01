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

using std::cout; using std::endl;

const double FMGEV = 5.08;

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
        double maxr = std::sqrt(2.0*B_p)*10; // Up to 10 Gaussian widths away, note that large r is very unlikely
        
        do{
            radius = gsl_rng_uniform(global_rng) * maxr;
        } while (gsl_rng_uniform(global_rng) > RadiusDistribution(radius));
        
        // Sample angle
        double costheta = 2.0*(gsl_rng_uniform(global_rng)-0.5);
        double sintheta = 2.0*(gsl_rng_uniform(global_rng)-0.5);     // TODO is this really uniform in theta
        Vec tmpvec(radius*costheta, radius*sintheta);
        quarks.push_back(tmpvec);
        quark_bp.push_back(B_q);  
    }
    
    
}

Ipsat_Proton::Ipsat_Proton()
{
    saturation=true;
    B_p = 4.0; // GeV^-2
    B_q = B_p/3.0;
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
    
    if (saturation)
        return 1.0 - std::exp( - r.LenSqr() * 1.0/quarks.size() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  );
    else
        return r.LenSqr() * 1.0/quarks.size() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  ;

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
double Ipsat_Proton::RadiusDistribution(double r)
{
    
    // We want to have on averate the same distribution as in the IPsat model
    return 0.5*std::exp( - r*r / (2.0*B_p));
    // Factor 0.5 should not matter, in just puts the probability always < 1 to accept
    // a sampled radius. Any factor ]0,1] should work (TEST!)
    
    /*
     * Exponential distribution, doesnt seemt to work with HERA
    r = r / 5.08;   // r to fm
    
    //a = \sqrt{12}/R_p = 3.87, with R_p = 0.895 from
    //http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1
    double a = std::sqrt(12)/R_p;
    double norm = 1.0 / ( std::pow(2.0/a, 2)*std::exp(-2) );    // Normalized to unity at maximum
    
    return norm * r*r*std::exp(-a*r);
     */
    
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

double Ipsat_Proton::Amplitude(double xpom, Vec q1, Vec q2)
{
    double quark[2] = {q1.GetX(), q1.GetY() };
    double antiquark[2] = {q2.GetX(), q2.GetY() };
    return Amplitude(xpom, quark, antiquark);
}