/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2022
 */

#ifndef diffraction_hpp
#define diffraction_hpp

#include "dipole.hpp"
#include "qcd.hpp"
#include "wave_function.hpp"


class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    double ScatteringAmplitude(double xpom, double Qsqr, double t,Polarization pol=T, bool real_part=true);
    double ScatteringAmplitudeF(double xpom, double Qsqr, double b, Polarization pol=T, bool real_part=true);
    double ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol=T, bool real_part=true);
    double ScatteringAmplitudeIntegrand2(
    double xpom, double Qsqr, double r, double theta_r,
    double b, double theta_b, double z, Polarization pol, bool real_part); 

    // Calculate scattering amplitude in case of cylinderical cymmetry (e.g. ipsat with no constituent quarks)
    double ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t, Polarization pol=T);
    double ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z, Polarization pol=T);
    
    
    double LogDerivative(double xpom, double Qsqr, double t, Polarization pol=T);   // der ln A / der y
    double Correction(double xpom, double Qsqr, double t, Polarization pol=T);
    
    void SetNumOfAverages(int n);
    
    void SetZLimit(double zl) { zlimit = zl; }

    bool ShowVegasIterations() { return show_vegas_iterations; }
    void ShowVegasIterations(bool s) { show_vegas_iterations = s; }
    
    DipoleAmplitude* GetDipole();
    WaveFunction* GetWaveFunction();
    
    // Set maximum dipole size
    void SetMaxR(double maxr) { MAXR = maxr; }
    double MaxR() { return MAXR; }
	

private:
    double MAXR;
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    int num_of_averages;
    
    double zlimit;
    
    double beta;
    double show_vegas_iterations;
    
};
#endif /* diffraction_hpp */
