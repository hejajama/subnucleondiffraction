/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a nucleus that consist of nucleons described ty IPsat
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "ipsat_nucleons.hpp"
#include "vector.hpp"
#include <gsl/gsl_rng.h>
#include <tools/config.hpp>
#include <sstream>

using namespace Amplitude;


double Ipsat_Nucleons::Amplitude(double xpom, double q1[2], double q2[2] )
{
    // Nucleon transveser coodinates are now nucleons[i].GetX() and GetY()
    Vec q(q1[0], q1[1]);
    Vec qbar(q2[0],q2[1]);
    
    Vec b = q + qbar;
    b = b*0.5;
    
    Vec r = q - qbar;
    
    // Need to calculate sum_i T_p(|b-b_i|), where b is the center of the dipole
    // and b_i is the center of the nucleon i
    // Note that at this point we do not take into account z, so b in geometric average, not
    // center of mass. Should we change this????
    
    // Currently should correspond to calculation arXiV:1011.1988
    
    double tpsum = 0;
    for (unsigned int i=0; i < nucleons.size(); i++)
    {
        Vec projection (nucleons[i].GetX(), nucleons[i].GetY());
        Vec deltab = b - projection;
        tpsum = tpsum + Tp(deltab.Len());
    }
    
    if (saturation)
        return 1.0 - std::exp( - r.LenSqr() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  );
    else
        return r.LenSqr() * gdist.Gluedist(xpom, r.LenSqr()) * tpsum  ;
}



/*
 * Initialize target nucelus
 * In this case: sample nucleon positions
 * Does this really sample from the WS distribution???
 */
void Ipsat_Nucleons::InitializeTarget()
{
    double smallestdist=999; double largestdist=0;
    double tmpdist;
    
    // Randomness
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    
    
    do
    {
        nucleons.clear(); smallestdist=999; largestdist=0;
        for (int i=0; i<A; i++)
        {
            Vec tmp;
            do {
                Vec tmpvec (2.0*(gsl_rng_uniform(r)-0.5)*MaxR(),
                            2.0*(gsl_rng_uniform(r)-0.5)*MaxR(),
                            2.0*(gsl_rng_uniform(r)-0.5)*MaxR());
                tmp=tmpvec;
            } while (gsl_rng_uniform(r) > WS_unnorm(tmp.Len())); // WS distribution!
            nucleons.push_back(tmp);
        }
        
        
        // Check smallestdist and largest dist
        for (int i=0; i<A; i++)
            for (int j=0; j<i; j++)
            {
                tmpdist = sqrt( SQR(nucleons[j].GetX() - nucleons[i].GetX())
                               + SQR(nucleons[j].GetY() - nucleons[i].GetY())
                               + SQR(nucleons[j].GetZ()-nucleons[i].GetZ()));
                if (tmpdist < smallestdist) smallestdist=tmpdist;
                if (tmpdist > largestdist) largestdist=tmpdist;
            }
        
        // Try again if smallestdist/largestdist does not match the limits
        if (smallestdist < mindist )
            cerr << "Again... Smallest dist: " << smallestdist << " , largest: " << largestdist << endl;
    } while (smallestdist < mindist );
    
    
    gsl_rng_free(r);
}


Ipsat_Nucleons::Ipsat_Nucleons()
{
    A=197;
    mindist = 0;
    ws_delta=0.54*FMGEV;
    ws_ra = 1.12 * std::pow(A, 1.0/3.0) * FMGEV;
    
    B_p = 4.0; // Teels the proton transverse size, GeV^2
    
    saturation = true;
    
}

// Max radius used to sample nucleons, depends on A, in GeV^-1
double Ipsat_Nucleons::MaxR()
{
    return 80;
}


double Ipsat_Nucleons::WS_unnorm(double r
                                 )
{
    return 1.0 / (1+exp((r-ws_ra)/ws_delta));
    
}

double Ipsat_Nucleons::Tp(double b)
{
    return 1.0/(2.0*M_PI*B_p)*std::exp(- b*b / (2.0*B_p));
}

std::vector<Vec> &Ipsat_Nucleons::GetNucleons()
{
    return nucleons;
}

void Ipsat_Nucleons::SetSaturation(double s)
{
    saturation = s;
}


std::string Ipsat_Nucleons::InfoStr()
{
    std::stringstream ss;
    ss << "IPsat model ";
    if (saturation)
        ss << "with saturation";
    else
        ss << "no saturation";
    
    ss <<". A=" << A << ", B_p=" << B_p << ", minimum nucleon-nucleon distance=" << mindist;
    return ss.str();
}