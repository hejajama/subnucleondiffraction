/*
 * Nucleus that consists of nucleons
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "nucleons.hpp"
#include "vector.hpp"
#include "ipsat_proton.hpp"
#include "subnucleon_config.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>

const double FMGEV = 5.08;

using std::cout;
using std::endl;

double Nucleons::Amplitude(double xpom, double q1[2], double q2[2] )
{
    // Calculate scattering amplitude in Glauber approach
    // Idea: the S matrix is a product of dipole-nucleon S matrixes
    // S(r,b) = prod_i S_p(r, b-b_i);
    // where b_i is the center of nucleon i
    
    double smat=1;
    Vec qv1(q1[0],q1[1]);
    Vec qv2(q2[0],q2[1]);
    for (int i=0; i<A; i++)
    {
        // Calculate the quark and antiquark coordinates in the frame where the nucleon i is at the origin
        Vec new_q1  = qv1 - nucleon_positions[i];
        Vec new_q2 = qv2 - nucleon_positions[i];
        smat = smat * (1.0 - nucleons[i].Amplitude(xpom, new_q1, new_q2));
        
    }
    return 1.0-smat;
}

void Nucleons::SetQuarkWidth(double bq)
{
    B_q = bq;
}

void Nucleons::SetProtonWidth(double bp)
{
    B_p = bp;
}

void Nucleons::InitializeTarget()
{
    nucleons.clear();
    nucleon_positions.clear();
    
    
    for (int i=0; i<A; i++)
    {   
    
        double maxr = 3.0*ws_ra;
        
        Vec tmp;
        do {
            Vec tmpvec (2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr);
            tmp=tmpvec;
        } while (gsl_rng_uniform(global_rng) > WS_unnorm(tmp.Len())); // WS distribution!
        nucleon_positions.push_back(tmp);

        Ipsat_Proton p(&gdist);
        p.SetQuarkWidth(B_q);
        p.SetProtonWidth(B_p);
        p.InitializeTarget();
        nucleons.push_back(p);
        

    }
    
}

Nucleons::Nucleons()
{
    A=197;
    ws_delta=0.54*FMGEV;
    ws_ra = 1.12 * std::pow(A, 1.0/3.0) * FMGEV;

}


double Nucleons::WS_unnorm(double r )
{
    return 1.0 / (1+exp((r-ws_ra)/ws_delta));
    
}