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

class DipoleAmplitude {
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    virtual double Amplitude(double xpom, double q1[2], double q2[2] );
    
    // Setup the target. In practice sample nucleon positions / quark positions / whatver.
    // This is called before Amplitude() is evaluated. When averaging over nuclear configurations
    // this always samples a new configuration.
    virtual void InitializeTarget();
    
    // Info string about the dipole
    virtual std::string InfoStr();
};

#endif