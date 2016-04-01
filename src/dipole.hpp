/*
 * Dipole cross section with explicit impact parameter and r dependence
 * Abstract class
 * All nuclear (target) dependence is described by this class]
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef dipole_h_
#define dipole_h_

#include <string>
#include "vector.hpp"

class DipoleAmplitude {
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    virtual double Amplitude(double xpom, double q1[2], double q2[2] );
    
    virtual double Amplitude(double xpom, Vec q1, Vec q2);
    
    virtual double AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] );
    
    double AmplitudeImaginaryPart(double xpom, Vec q1, Vec q2 );
    
    virtual void SetSkewedness(bool s); // Enalbe/disable skewedness, only used by ipsat
    
    // Setup the target. In practice sample nucleon positions / quark positions / whatver.
    // This is called before Amplitude() is evaluated. When averaging over nuclear configurations
    // this always samples a new configuration.
    virtual void InitializeTarget();
    
    // Info string about the dipole
    virtual std::string InfoStr();
    
    // Calculate saturation scale
    // Note that in this general case we actually do not have an unique way to define Q_s
    // Here we define the saturation scale around point b via r_s such that
    // N(b + 0.5 r_s, b - 0.5 r_s) = 1 - exp(-0.5), and
    // r_s^2 = 2/Q_s^2
    // That is, we choose the arbitrary direction to be along the x axis
    double SaturationScale(double xpom, Vec b);
    
    virtual double Density(Vec b);  // T_p(b), or whatever replaces T_b
};

#endif