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
        tpsum = tpsum + Tp(deltab.Len(), i);
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
    
    cout << "# Sampling IPsat nucleus, rng name " << gsl_rng_name(r) << " seed " << gsl_rng_default_seed << " first value " << gsl_rng_get(r) << " first float " << gsl_rng_uniform(r) << endl;
    
    
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
    
    
    // Sample size of the nucleon
    if (radius_fluctuation_fraction > 0.0001)
    {
        for (int i=0; i<A; i++)
        {
            double r0 = std::sqrt(2.0*B_p);
            
            // How small/large distances can in principle be sampled, max deltar
            // is r0*maxradius_n. However, the large/small radiuses are extremely unlikely
            double maxradius_n = 5;
            
            double deltar = r0*radius_fluctuation_fraction;
            double new_r=0;
            
            do
            {
                double rand = 2.0*(gsl_rng_uniform(r)-0.5);
                new_r = r0 + rand*r0*maxradius_n;
            } while (new_r <= 0 or gsl_rng_uniform(r) > NucleonRadiusDistribution(r0, new_r, deltar));
            
            B_ps[i] = new_r * new_r / 2.0;
            
        }
            
    }
    

    
    gsl_rng_free(r);
}


Ipsat_Nucleons::Ipsat_Nucleons()
{
    A=197;
    mindist = 0;
    ws_delta=0.54*FMGEV;
    ws_ra = 1.12 * std::pow(A, 1.0/3.0) * FMGEV;
    
    B_p = 4.0; // Tells the proton transverse size, GeV^-2
    
    for (int i=0; i<A; i++)
    {
        // If fluctuating nucleon sizes are used, this is overwritten at InitializeTarget()
        B_ps.push_back(B_p);
        
    }
    radius_fluctuation_fraction = 0;
    
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

// Proton thicness, by defualt proton_id is -1, in that case use the default B_p
double Ipsat_Nucleons::Tp(double b, int proton_id)
{
    double bp = B_p;
    if (proton_id >=0 )
        bp = B_ps[proton_id];
    return 1.0/(2.0*M_PI*bp)*std::exp(- b*b / (2.0*bp));
}

std::vector<Vec> &Ipsat_Nucleons::GetNucleons()
{
    return nucleons;
}

std::vector<double> &Ipsat_Nucleons::GetB_ps()
{
    return B_ps;
}

void Ipsat_Nucleons::SetSaturation(double s)
{
    saturation = s;
}

void Ipsat_Nucleons::SetFluctuatingNucleonSize(double f)
{
    radius_fluctuation_fraction = f;
}

double Ipsat_Nucleons::NucleonRadiusDistribution(double r0, double r, double width)
{
    return std::exp(- SQR(r-r0)/(2.0*SQR(width)));
};

std::string Ipsat_Nucleons::InfoStr()
{
    std::stringstream ss;
    ss << "IPsat model ";
    if (saturation)
        ss << "with saturation";
    else
        ss << "no saturation";
    
    ss <<". A=" << A << ", B_p=" << B_p << ", minimum nucleon-nucleon distance=" << mindist;
    ss << ". Fluctuation parameter radius_fluctuation_fraction: " << radius_fluctuation_fraction;
    return ss.str();
}