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
#include <complex>

class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    double ScatteringAmplitude(double xpom, double Qsqr, double t,Polarization pol=T );
    double ScatteringAmplitude_at_fixed_b(double xpom, double Qsqr, double b, double theta_b, Polarization pol);
    double ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol=T);
    double ScatteringAmplitudeIntegrand_fixed_b(double xpom, double Qsqr,  double r, double theta_r, double b, double theta_b, double z, Polarization pol);
    double Relative_P_abd_B_Inte(double mv, double root_snn, double theta_BigP, int Z, double t, double RA, double Low, double High, bool DacayToScalar);
    double Relative_P_abd_B_Inte_mc(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, double theta_BigP, int Z, 
                                    double t, double RA, double Low, double High, bool DacayToScalar);

    double b_Inte(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, double B, double theta_B, 
                  int Z, double t, double BigP, double theta_BigP, int M12reim);
                           
    double Soft_photon_ScatteringAmplitude(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, 
           double theta_BigP, int Z, double t, double RA, double Low, double High, double daughter_mass, bool DacayToScalar);
    std::pair<double, double> ScatteringAmplitude_noexp_Integrand(double xpom, double Qsqr, double t, double r, double theta_r, 
                              double b, double theta_b, double bp, double theta_bp, double bp2, double theta_bp2, double z, Polarization pol);

    // Calculate scattering amplitude in case of cylinderical cymmetry (e.g. ipsat with no constituent quarks)
    double ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t, Polarization pol=T);
    double ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z, Polarization pol=T);
    
    std::complex<double> ScatteringAmplitudeIntegrand_reim(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol);
    double LogDerivative(double xpom, double Qsqr, double t, Polarization pol=T);   // der ln A / der y
    double Correction(double xpom, double Qsqr, double t, Polarization pol=T);
    
    void SetNumOfAverages(int n);
    
    void SetZLimit(double zl) { zlimit = zl; }
    
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
    
};
#endif /* diffraction_hpp */
