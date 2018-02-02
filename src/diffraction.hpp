/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef diffraction_hpp
#define diffraction_hpp

#include "dipole.hpp"
#include "wave_function.hpp"
#include <tools/config.hpp>

using namespace Amplitude;

enum DIJET_COMPONENT
{
    X,
    Y
};

class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    double ScatteringAmplitude(double xpom, double Qsqr, double t,Polarization pol=T );
    double ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol=T);
    
    // Calculate scattering amplitude in case of cylinderical cymmetry (e.g. ipsat with no constituent quarks)
    double ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t, Polarization pol=T);
    double ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z, Polarization pol=T);
    
    
    double LogDerivative(double xpom, double Qsqr, double t, Polarization pol=T);   // der ln A / der y
    double Correction(double xpom, double Qsqr, double t, Polarization pol=T);
    
    void SetNumOfAverages(int n);
    
    void SetZLimit(double zl) { zlimit = zl; }
    
    DipoleAmplitude* GetDipole();
    WaveFunction* GetWaveFunction();
    
    void SetDijetComponent(DIJET_COMPONENT c) { dijet_component=c; }
    DIJET_COMPONENT GetDijetComponent() { return dijet_component; }
private:
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    int num_of_averages;
    
    double zlimit;
    
    double beta;
    
    DIJET_COMPONENT dijet_component;
    
};
#endif /* diffraction_hpp */
