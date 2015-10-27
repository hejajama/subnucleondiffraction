/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for dipole-smooth nucleus scattering for testing/comparisons
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "smooth_ws_nuke.hpp"
#include <tools/tools.hpp>
#include <cmath>

using Amplitude::SQR;



using namespace std;

Smooth_ws_nuke::Smooth_ws_nuke()
{
    A=197;
    InitializeWSDistribution(A);
    
}

double Smooth_ws_nuke::Amplitude(double xpom, double q1[2], double q2[2] )
{
        
    // Very crude toy model: neglet different positions for quarks, only take impact parameter dependence from the Woods Saxon
    // This is very roughly like the ipsat model where we replace proton T_p by T_A, and drop all the constants at this point
    // Note: close to center this is roughly one and drops to zero at the edges
    
    double r = std::sqrt( SQR(q1[0]-q2[0]) + SQR(q1[1]-q2[1]) );
    
    // Take the nucleon density at the geometric center of the two quarks
    double b = std::sqrt( SQR( (q1[0]+q2[0])/2.0 ) + SQR( (q1[1]+q2[1])/2.0) );
    
    return 1.0 - exp( -A*gdist.Gluedist(xpom, r*r) *r*r* T_A(b, A));
    
    // As T_A is normalized to unity, we get no extra normalizatino factor for it
    
}
