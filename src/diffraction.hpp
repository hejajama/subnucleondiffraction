/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef diffraction_hpp
#define diffraction_hpp

#include "dipole.hpp"
#include "wave_function.hpp"

enum Polarization
{
    TRANSVERSE,
    LONGITUDINAL
};

class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    double ScatteringAmplitude(double xpom, double Qsqr, double t,Polarization pol=TRANSVERSE );
    double ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol=TRANSVERSE);
    
    // Calculate scattering amplitude in case of cylinderical cymmetry (e.g. ipsat with no constituent quarks)
    double ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t);
    double ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z);
    
    
    double LogDerivative(double xpom, double Qsqr, double t);   // der ln A / der y
    double Correction(double xpom, double Qsqr, double t);
    
    void SetNumOfAverages(int n);
    
    DipoleAmplitude* GetDipole();
    WaveFunction* GetWaveFunction();
private:
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    int num_of_averages;
    
    double beta;
    
};
#endif /* diffraction_hpp */
