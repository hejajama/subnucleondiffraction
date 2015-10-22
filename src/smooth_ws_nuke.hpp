/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for dipole-smooth nucleus scattering for testing/comparisons
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef smooth_ws_nuke_hpp
#define smooth_ws_nuke_hpp

#include "dipole.hpp"

class Smooth_ws_nuke : public DipoleAmplitude
{
public:
    Smooth_ws_nuke();
    double Amplitude(double xpom, double q1[2], double q2[2] );
private:
    int A;
    
};

#endif /* smooth_ws_nuke_hpp */
