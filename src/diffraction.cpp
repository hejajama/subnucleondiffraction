/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#include "diffraction.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
}

/* 
 * Diffractive scattering amplitude
 * t: squared momentum transfer
 * xpom: Bjorken x
 * Qsqr: Q^2
 *
 * Calculated by integrating over the transverse positions of b and r (2d vecs, use monte carlo) and momentum fraction z
 */

struct Inthelper_amplitude
{
    Diffraction* diffraction;
    double xpom;
    double Qsqr;
    double t;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

double Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;

    // Do MC integral over impact parameters and dipole sizes
    
    // Currently hardcoded parameters for jpsi and gold:
    // Impact parameter up to 100 GeV^-1
    // Dipole size up to 10 GeV^-1
    // MC integration parameters: b, theta_b, r, theta_r
    double lower[4] = {0,0,0,0};
    double upper[4] = {100, 2.0*M_PI, 10, 2.0*M_PI};
    
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 4;
    F.params = &helper;
    
}


double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par)
{
    
    return 0;
}

// Recall quark and gluon positions:
// Quark: b + zr
// Antiquark: b - (1-z) r

// If do like Lappi, Mantysaari: set z=1/2 here
/*
double bx = b*cos(theta_b);
double by = b*sin(theta_b);
double rx = r*cos(theta_r);
double ry = r*sin(theta_r);

// q and antiq positions
double qx = bx + z*rx; double qy = by + z*ry;
double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
*/


// Helpers
DipoleAmplitude
* Diffraction::GetDipole()
{
    return dipole;
}
WaveFunction* Diffraction::GetWaveFunction(){
    return wavef;
}