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
#include <vector>
#include <complex>


class Diffraction
{
public:
    Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_);
    
    // Calculate amplitude A, this will later be averaged and squared
    std::complex<double> ScatteringAmplitude(double xpom, double Qsqr, double t,
        Polarization pol=T);
    // Forward amplitude at t=0 integrating over theta_b internally (Suave vector integration)
    std::complex<double> ScatteringAmplitudeF(double xpom, double Qsqr, double b,
        double theta_b, Polarization pol=T, double* integrand_mod_sqr=nullptr);

    std::complex<double> ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t,
        double r, double theta_r, double b, double theta_b, double z,
        Polarization pol=T);

    // Calculate scattering amplitude in case of cylinderical cymmetry (e.g. ipsat with no constituent quarks)
    double ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, 
        double t, Polarization pol=T);
    double ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, 
        double Qsqr, double t, double r, double b, double z, Polarization pol=T);
    
    
    double LogDerivative(double xpom, double Qsqr, double t, Polarization pol=T);   // der ln A / der y
    double Correction(double xpom, double Qsqr, double t, Polarization pol=T);
       
    void SetZLimit(double zl) { zlimit = zl; }

    DipoleAmplitude* GetDipole();
    WaveFunction* GetWaveFunction();
    
    // Set maximum dipole size
    void SetMaxR(double maxr) { MAXR = maxr; }
    double MaxR() { return MAXR; }

    struct TotalCrossSectionData {
        std::vector<double> b;      // GeV^-1 centers, size nbperp
        std::vector<double> theta;  // angles, size ntheta
        // The following vectors are stored on a flattened (b,theta) grid
        // with index i = ib*ntheta + itheta.
        std::vector<std::complex<double>> F_T;            // size nbperp*ntheta
        std::vector<std::complex<double>> F_L;            // empty if Q^2=0
        std::vector<double> F_T_sqr;                      // |F_T|^2 per (b,theta)
        std::vector<double> F_L_sqr;                      // |F_L|^2 per (b,theta), empty if Q^2=0
        // Integral over internal variables of |integrand|^2 at fixed (b,theta)
        std::vector<double> F_T_integrand_sqr;            // size nbperp*ntheta
        std::vector<double> F_L_integrand_sqr;            // empty if Q^2=0
    };

    TotalCrossSectionData ComputeTotalCrossSection(double xpom, double Qsqr,
        int nbperp, double maxb, int ntheta);
	

private:
    double MAXR;
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    
    double zlimit;
    
    double beta;
    double show_vegas_iterations;
    
};
#endif /* diffraction_hpp */
