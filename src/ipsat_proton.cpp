/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "ipsat_proton.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>

using std::cout; using std::endl;

void Ipsat_Proton::InitializeTarget()
{
    double smallestdist=999; double largestdist=0;
    double tmpdist;
    
    // Randomness
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    
    cout << "# Sampling IPsat proton, rng name " << gsl_rng_name(r) << " seed " << gsl_rng_default_seed << " first value "; cout << gsl_rng_get(r); cout << " first float " << gsl_rng_uniform(r) << endl;
    
    
    quarks.clear();
    quark_radii.clear();
    
    // Sample 3 quarks
    for (int i=0; i<3; i++)
    {
        // Radius from uniform distribution
        double radius=0;
        
        do{
            radius = gsl_rng_uniform(r) * maxr;
        } while (gsl_rng_uniform(r) > RadiusDistribution(radius));
        
        // Sample angle
        double costheta = 2.0*(gsl_rng_uniform(r)-0.5);
        double sintheta = 2.0*(gsl_rng_uniform(r)-0.5);     // TODO is this really uniform in theta
        Vec tmpvec(radius*costheta, radius*sintheta);
        quarks.push_back(tmpvec);
        quark_radii.push_back(R_p / 3.0);   // Quark radius is 1/3 proton radius
    }
    
    
}

Ipsat_Proton::Ipsat_Proton()
{
    saturation=true;
    R_p = 0.895; // 0.895 fm is the proton radius
    
    maxr = 5.0 * R_p * 5.08; // Max radius sampled in GeV^-1
    // factor 5 in principle allows very large protons, but they are very unlikely
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
        return 1.0 - std::exp( - r.LenSqr() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  );
    else
        return r.LenSqr() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  ;

}

std::vector<Vec> &Ipsat_Proton::GetQuarks()
{
    return quarks;
}

std::vector<double> &Ipsat_Proton::GetRadii()
{
    return quark_radii;
}

double Ipsat_Proton::QuarkThickness(double r, int i)
{
    double bp = quark_radii[i]*quark_radii[i]/2.0;
    return 1.0/(2.0*M_PI*bp)*std::exp(- r*r / (2.0*bp));
}

/*
 * Quark distances from the origin are sampled from this distribution
 */
double Ipsat_Proton::RadiusDistribution(double r)
{
    r = r / 5.08;   // r to fm
    
    //a = \sqrt{12}/R_p = 3.87, with R_p = 0.895 from
    //http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1
    double a = std::sqrt(12)/R_p;
    double norm = 1.0 / ( std::pow(2.0/a, 2)*std::exp(-2) );    // Normalized to unity at maximum
    
    return norm * r*r*std::exp(-a*r);
    
}