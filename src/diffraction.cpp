/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#include "diffraction.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include "subnucleon_config.hpp"



Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    num_of_averages = 1;
}

/*
 * Calculate total diffractive cross section
 * Requires average over different nucleon configurations
 *
 * First take the square, then average
 */
double Diffraction::TotalCrossSection(double xpom, double Qsqr, double t)
{
    // Testing: only amplitude squared
    double amp = ScatteringAmplitude(xpom, Qsqr, t);
    return amp*amp/(16.0*M_PI);
}

/*
 * Calculate total diffractive cross section
 * Requires average over different nucleon configurations
 *
 * First average, then take the square
 */
double Diffraction::CoherentCrossSection(double xpom, double Qsqr, double t)
{
    
    return 0;
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
    double r;
    double theta_r;
    double b;
    double theta_b;
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
    
    double result,error;
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_monte_miser_state *s = gsl_monte_miser_alloc(4);
    gsl_monte_miser_integrate(&F, lower, upper, 4, MCINTPOINTS, r, s, &result, &error);
    gsl_monte_miser_free(s);
    gsl_rng_free(r);
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    return result;
    
}



double Inthelperf_amplitude_z(double z, void* p);

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->b = vec[0];
    helper->theta_b=vec[1];
    helper->r = vec[2];
    helper->theta_r = vec[3];
    
    // Use the following if we integrate over z using gsl_integration
    // For now we neglect (1-z)rDelta term, so can only integrate the wave function overlap
    // over z. Thus, at this point we just call the integrand and in the integrand integrate the overlap over z
    // Note that if integrand is changed, also this must be changed
    
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, 0.5); // Put z=0.5 as it sets b to the geometric average of quarks
    /*
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(ZINT_INTERVALS);
    
    gsl_function F;
    F.function =
    &Inthelperf_amplitude_z;
    F.params = par;
    double result,error;
    int status = gsl_integration_qags(&F, 0, 1, 0, ZINT_RELACCURACY, ZINT_INTERVALS, w, &result, &error);
    
    if (status)
        cerr << "#ZINT failed, result " << result << " relerror " << error << " r " << helper->r << " b " << helper->b << endl;
    
    gsl_integration_workspace_free(w);
    
    return result;
     */
}

double Inthelperf_amplitude_z(double z, void* p)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)p;
    
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z);
}

double Diffraction::ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z)
{ 
    
    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    
    // If do like Lappi, Mantysaari: set z=1/2 here
    
    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    
    // q and antiq positions
    double tmpz = z;
    z=0.5;      // Currently use b as geometric average
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    z=tmpz;
    
    double res = 0;
    
    res = 2.0 * r * b;
    
    //res *= wavef->PsiSqr_tot(Qsqr, r, z)/(4.0*M_PI); // Wavef
    // As this integrand is now not integrated over z
    // Also, this only works for photoproduction!
    res *= wavef->PsiSqr_T_intz(Qsqr, r) / (4.0*M_PI);
    if (Qsqr>0)
        cerr << "ScatteringAmplitudeIntegrand currently only works for photoproduction!" << endl;
    
    double delta = std::sqrt(t);
    
    // Real part cos, imaginary part -sin
    //res *= std::cos( b*delta*std::cos(theta_b) - (1.0 - z)*r*delta*std::cos(theta_r));
    res *= std::cos( b*delta*std::cos(theta_b));    // Neglect z now
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    
    res *= dipole->Amplitude(xpom, x1, x2 );
    
    return res;
    

}



// Helpers

void Diffraction::SetNumOfAverages(int n)
{
    num_of_averages = n;
}

DipoleAmplitude
* Diffraction::GetDipole()
{
    return dipole;
}
WaveFunction* Diffraction::GetWaveFunction(){
    return wavef;
}