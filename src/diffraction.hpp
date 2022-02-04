/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef diffraction_hpp
#define diffraction_hpp

#include "dipole.hpp"
#include <amplitudelib/wave_function.hpp>
#include <tools/interpolation.hpp>
#include <tools/config.hpp>
#include <string>
#include <vector>

enum COMPONENT
{
    X,
    Y,
    M0,
    M1x,
    M1y
};

using namespace Amplitude;

class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    ~Diffraction();
    
    // Calculate amplitude A, this will later be averaged and squared
    std::vector<double> ScatteringAmplitude(double xpom, double Qsqr, double t, double B, double theta_B, bool real_part, std::vector<COMPONENT> complist, Polarization pol=T);
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
    
    // Set maximum dipole size
    void SetMaxR(double maxr) { MAXR = maxr; }
    double MaxR() { return MAXR; }
    
    void InitializeIsInterpolator(std::string datafile);
    Interpolator *GetIsInterpolator() { return Is_interpolator; }
    

	

private:
    double MAXR;
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    int num_of_averages;
    
    double zlimit;
    
    double beta;

    Interpolator *Is_interpolator; // Function I_s(B) needed with finite photon kT
    
};
#endif /* diffraction_hpp */
