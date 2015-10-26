/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef diffraction_hpp
#define diffraction_hpp

#include "dipole.hpp"
#include "wave_function.hpp"

class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    double ScatteringAmplitude(double xpom, double Qsqr, double t);
    double ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z);
    
    double CoherentCrossSection(double xpom, double Qsqr, double t);
    double TotalCrossSection(double xpom, double Qsqr, double t);
    
    void SetNumOfAverages(int n);
    
    DipoleAmplitude* GetDipole();
    WaveFunction* GetWaveFunction();
private:
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    int num_of_averages;
    
};
#endif /* diffraction_hpp */
