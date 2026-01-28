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
    
    // Scattering amplitude, integrated over b, r, z
    // See docs/t_integrated_vm_production for definition
    // Note: t is taken to be positive
    std::complex<double> ScatteringAmplitude(double xpom, double Qsqr, double t,
        Polarization pol=T);
    std::complex<double> ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t,
        double r, double theta_r, double b, double theta_b, double z,
        Polarization pol=T);    

    // Scattering amplitude integrated over t at fixed b
    // See docs/t_integrated_vm_production for definition
    std::complex<double> ScatteringAmplitude_tIntegrated(double xpom, double Qsqr, double b,
        double theta_b, Polarization pol=T);

    

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
    };

    TotalCrossSectionData ComputeTotalCrossSection(double xpom, double Qsqr,
        int nbperp, double maxb, int ntheta);
    
    void SetFactorizeZInt(bool fz) { factorize_zint = fz; } 
    bool GetFactorizeZInt() const { return factorize_zint; }
    void SetMCIntPoints(unsigned int points) { mcintpoints = points; }
    unsigned int GetMCIntPoints() const { return mcintpoints; }

private:
    double MAXR;
    DipoleAmplitude* dipole;
    WaveFunction* wavef;
    
    double zlimit;
    
    double beta;

    bool factorize_zint = false;

    unsigned int mcintpoints = 1e5;
    
};
#endif /* diffraction_hpp */
