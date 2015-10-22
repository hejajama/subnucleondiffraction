/*
 * Dipole cross section with explicit impact parameter and r dependence
 * Abstract class
 * All nuclear (target) dependence is described by this class]
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef dipole_h_
#define dipole_h_

class DipoleAmplitude {
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    virtual double Amplitude(double xpom, double q1[2], double q2[2] );
};

#endif