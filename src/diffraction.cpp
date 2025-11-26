/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015-2025
 */
#include "diffraction.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "subnucleon_config.hpp"
#include "nrqcd_wf.hpp"

using namespace std;

#include <complex>

#include <functional>
#include <iomanip>
#include <cstdlib>

#include "Cuba-4.2.2/cuba.h"


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    zlimit=0.00000001;
	MAXR=10*5.068;
    show_vegas_iterations=true;
}




/*
 * Derivative of the ampltiude, only for rotationally symmetric amplitude
 * Calculate d ln N / d y, y = ln 1/x
 */
struct AmplitudeDerHeler
{
    Diffraction* diff;
    Polarization pol;
    double Qsqr;
    double t;
};

double AmplitudeDerHelperf(double y, void* p)
{
    AmplitudeDerHeler* par = (AmplitudeDerHeler*)p;
    double x = exp(-y);
    double res = std::log(std::abs(par->diff->ScatteringAmplitudeRotationalSymmetry(x, par->Qsqr, par->t, par->pol)));
    return res;
}
double Diffraction::LogDerivative(double xpom, double Qsqr, double t, Polarization pol)
{
    gsl_function F;
    F.function=&AmplitudeDerHelperf;
    AmplitudeDerHeler par; par.Qsqr=Qsqr; par.t=t; par.pol=pol;
    par.diff  = this;
    F.params = &par;
    double result,abserr;
    double y = std::log(1.0/xpom);
    gsl_deriv_central (&F, y, 0.1 , &result, &abserr);
    
    //cout << "Der " << result << " err " << abserr << endl;
    return result;
}

/* Calculate total correction
 */
double Diffraction::Correction(double xpom, double Qsqr, double t, Polarization pol)
{
    double lambda = LogDerivative(xpom, Qsqr, t, pol);
    
    double beta = std::tan(lambda*M_PI/2.0);
    
    double Rg = std::pow(2.0, 2.0*lambda+3)/std::sqrt(M_PI) * gsl_sf_gamma(lambda+5.0/2.0)/gsl_sf_gamma(lambda+4.0);
    return (1.0+beta*beta)*Rg*Rg;
}


/* 
 * Diffractive scattering amplitude
 * t: squared momentum transfer
 * xpom: Bjorken x
 * Qsqr: Q^2
 *
 * Calculated by integrating over the transverse positions of b and r (2d vecs, use monte carlo) and momentum fraction z
 */
// Full t-dependent scattering amplitude using Suave over (b, r, theta_b, theta_r [, z])
std::complex<double> Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, Polarization pol) {
    struct SAParams {
        Diffraction* diff;
        double xp;
        double Q2;
        double t;
        double bmax;
        double rmin;
        double rmax;
        double zmin;
        bool fact;
        Polarization pol;
    };
    SAParams p{
        this,
        xpom,
        Qsqr, 
        t,
        10*5.068,
        1e-10,
        MAXR,
        zlimit,
        FACTORIZE_ZINT,
        pol
    };
    auto integrand = [](const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* ud)->int {
        SAParams* prm = static_cast<SAParams*>(ud);
        const double twoPi = 2.0*M_PI;
        const bool fact = prm->fact;
        // Unpack Cuba randoms
        const double xb = x[0];                // b in [0,1]
        const double xu = x[1];                // u=ln r in [0,1]
        const double xtb = x[2];               // theta_b in [0,1]
        const double xtr = x[3];               // theta_r in [0,1]
        const double xz = fact ? 0.5 : x[4];   // z in [0,1]
        // Map to physical
        const double b = prm->bmax * xb;
        const double umin = std::log(prm->rmin), umax = std::log(prm->rmax);
        const double u = umin + (umax-umin) * xu;
        const double r = std::exp(u);
        const double theta_b = twoPi * xtb;
        const double theta_r = twoPi * xtr;
        const double z = fact ? 0.5 : (prm->zmin + (1.0 - 2.0*prm->zmin) * xz);
        // Overlap and scalar prefactor (2 r b * overlap)
        double scalar = 2.0 * r * b;
        if (fact) {
            if (prm->diff->wavef->WaveFunctionType() == "NRQCD") {
                double delta_local = std::sqrt(prm->t);
                if (prm->pol == T)
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_T_intz(prm->Q2, r, delta_local, theta_r);
                else
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_L_intz(prm->Q2, r, delta_local, theta_r);
            } else {
                if (prm->pol == T)
                    scalar *= prm->diff->wavef->PsiSqr_T_intz(prm->Q2, r);
                else
                    scalar *= prm->diff->wavef->PsiSqr_L_intz(prm->Q2, r);
            }
        } else {
            const double inv4pi = 1.0/(4.0*M_PI);
            if (prm->pol == T)
                scalar *= prm->diff->wavef->PsiSqr_T(prm->Q2, r, z) * inv4pi;
            else
                scalar *= prm->diff->wavef->PsiSqr_L(prm->Q2, r, z) * inv4pi;
        }
        // Quark positions: factorized => z=1/2, non-factorized use (1-z) and z
        const double cos_tb = std::cos(theta_b);
        const double sin_tb = std::sin(theta_b);
        const double cos_tr = std::cos(theta_r);
        const double sin_tr = std::sin(theta_r);
        double bx = b * cos_tb;
        double by = b * sin_tb;
        double rx = r * cos_tr;
        double ry = r * sin_tr;
        double qx, qy, qbarx, qbary;
        if (fact) {
            qx = bx + 0.5*rx; qy = by + 0.5*ry;
            qbarx = bx - 0.5*rx; qbary = by - 0.5*ry;
        } else {
            qx = bx + (1.0 - z)*rx; qy = by + (1.0 - z)*ry;
            qbarx = bx - z*rx; qbary = by - z*ry;
        }
        double x1[2] = {qx,qy}; double x2[2] = {qbarx,qbary};
        std::complex<double> amp = prm->diff->dipole->ComplexAmplitude(prm->xp, x1, x2);
        // Phase factor with momentum transfer delta
        const double delta = std::sqrt(prm->t);
        if (delta > 0) {
            double phi;
            if (fact) {
                phi = b*delta*cos_tb;
            } else {
                phi = b*delta*cos_tb - (0.5 - z)*r*delta*cos_tr;
            }
            const std::complex<double> exponent = std::exp(std::complex<double>(0.0, -phi));
            amp *= exponent;
        }
        std::complex<double> val = scalar * amp;
        // Jacobian: bmax * (umax-umin) * (2pi)^2 * z-width (if unfactorized)
        const double J = prm->bmax * (umax-umin) * (twoPi*twoPi) * (fact ? 1.0 : (1.0 - 2.0*prm->zmin));
        const double measure_r = r; // from dr = r du
        f[0] = J * measure_r * static_cast<cubareal>(val.real());
        f[1] = J * measure_r * static_cast<cubareal>(val.imag());
        return 0;
    };
    const int ndim = FACTORIZE_ZINT ? 4 : 5;
    const int ncomp = 2;
    int nregions=0, neval=0, fail=0; double integral[2], error[2], prob[2];
    const int nvec = 1; const double epsrel = MCINTACCURACY, epsabs = 0.0;
    const int flags = 0, seed = 0;
    const int mineval = MCINTPOINTS/10; const int maxeval = MCINTPOINTS;
    const int nnew = mineval/20, nmin = 300; const double flatness = 1.0;
    Suave(ndim, ncomp, integrand, &p, nvec, epsrel, epsabs, flags, seed,
        mineval, maxeval, nnew, nmin, flatness,
        NULL, NULL, &nregions, &neval, &fail, integral, error, prob);
    return std::complex<double>(integral[0], integral[1]);
}

std::complex<double> Diffraction::ScatteringAmplitudeF(
    double xpom, double Qsqr, double b, Polarization pol, double* integrand_mod_sqr) {
    struct SuaveParams {
        Diffraction* diff;
        double xpom;
        double Q2;
        double b;
        double zmin;
        double rmin;
        double rmax;
        bool factorize;
        Polarization pol;
    };
    SuaveParams p{
        this,
        xpom,
        Qsqr,
        b,
        zlimit,
        1e-10,
        MAXR,
        FACTORIZE_ZINT,
        pol
    };
    auto integrand = [](const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *ud)->int{
        SuaveParams* prm = static_cast<SuaveParams*>(ud);
        const double twoPi = 2.0*M_PI;
        const bool fact = prm->factorize;
        const double umin = std::log(prm->rmin), umax = std::log(prm->rmax);
        const double xb = x[0];
        const double xr = x[1];
        const double xu = x[2];
        const double xz = fact ? 0.5 : x[3];
        const double theta_b_int = twoPi * xb;
        const double theta_r = twoPi * xr;
        const double u = umin + (umax-umin) * xu;
        const double r = std::exp(u);
        const double z = fact ? 0.5 : (prm->zmin + (1.0 - 2.0*prm->zmin) * xz);
        // Common factors
        // r from Jacobians, 2 as we have written sigma_qq = 2 N
        double scalar = 2.0 * r; // r from Jacobian (du->dr adds r)
        if (fact) {
            if (prm->diff->wavef->WaveFunctionType() == "NRQCD") {
                double delta = 0.0; // t=0
                if (prm->pol == T)
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_T_intz(prm->Q2, r, delta, theta_r);
                else
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_L_intz(prm->Q2, r, delta, theta_r);
            } else {
                if (prm->pol == T)
                    scalar *= prm->diff->wavef->PsiSqr_T_intz(prm->Q2, r);
                else
                    scalar *= prm->diff->wavef->PsiSqr_L_intz(prm->Q2, r);
            }
        } else {
            const double inv4pi = 1.0/(4.0*M_PI);
            if (prm->pol == T)
                scalar *= prm->diff->wavef->PsiSqr_T(prm->Q2, r, z) * inv4pi;
            else
                scalar *= prm->diff->wavef->PsiSqr_L(prm->Q2, r, z) * inv4pi;
        }
        // Recall quark and gluon positions:
        // Quark: b + zr
        // Antiquark: b - (1-z) r
        // If do like Lappi, Mantysaari: set z=1/2 here
        // Dipole amplitude
        const double cos_tb = std::cos(theta_b_int);
        const double sin_tb = std::sin(theta_b_int);
        const double cos_tr = std::cos(theta_r);
        const double sin_tr = std::sin(theta_r);
        double bx = prm->b * cos_tb;
        double by = prm->b * sin_tb;
        double rx = r * cos_tr;
        double ry = r * sin_tr;
        // q and antiq positions
        // Note my convention is that b is the center of the dipole (geometric center), not center of mass (z weighted)
        // Consequently I get (0.5-z)r. Delta phase
        double qx = bx + (1. - z) * rx;
        double qy = by + (1. - z) * ry;
        double qbarx = bx - z * rx;
        double qbary = by - z * ry;
        double x1[2] = {qx,qy};
        double x2[2] = {qbarx,qbary};
        std::complex<double> amp = prm->diff->dipole->ComplexAmplitude(prm->xpom, x1, x2);
        const double amp_r = amp.real();
        const double amp_i = amp.imag();
        // Overall Jacobian (theta_b, theta_r, u, z (if not factorized))
        const double J = (twoPi * twoPi) * (umax-umin) * (fact ? 1.0 : (1.0 - 2.0*prm->zmin));
        // Jacobian pieces:
        //  theta_b: 2pi, theta_r: 2pi  (in J)
        //  u = ln r mapping: u = umin + (umax-umin)*xu gives width (umax-umin) in J and dr = r du adds extra r
        //  optional z: width (1 - 2 zmin) in J
        // scalar currently includes factor 2*r * overlap; needs extra r from dr=r du
        const double measure_r = r; // from dr = r du
        // Components: real, imag, |integrand|^2 integral element (int |A|^2 dvars)
        f[0] = J * measure_r * scalar * amp_r; // real part
        f[1] = J * measure_r * scalar * amp_i; // imag part
        f[2] = J * measure_r * scalar * scalar * (amp_r*amp_r + amp_i*amp_i);
        return 0;
    };

    const int ndim = FACTORIZE_ZINT ? 3 : 4;
    const int ncomp = 3; // real, imag, |integrand|^2
    int nregions=0, neval=0, fail=0;
    double integral[3], error[3], prob[3];
    const int nvec = 1;
    const double epsrel = MCINTACCURACY, epsabs = 0.0;
    const int flags = 0, seed = 0;
    const int mineval = MCINTPOINTS/10; const int maxeval = MCINTPOINTS;
    const int nnew = mineval/20, nmin = 300; const double flatness = 1.0;
    Suave(ndim, ncomp, integrand, &p, nvec, epsrel, epsabs, flags, seed,
        mineval, maxeval, nnew, nmin, flatness,
        NULL, NULL, &nregions, &neval, &fail, integral, error, prob);
    if (integrand_mod_sqr)
        *integrand_mod_sqr = integral[2];
    return std::complex<double>(integral[0], integral[1]);
}

Diffraction::TotalCrossSectionData Diffraction::ComputeTotalCrossSection(
    double xpom, double Qsqr, int nbperp, double maxb) {
    TotalCrossSectionData out;
    out.b.resize(nbperp);
    out.F_T.assign(nbperp, std::complex<double>(0.,0.));
    if (Qsqr > 0) out.F_L.assign(nbperp, std::complex<double>(0.,0.));
    out.F_T_sqr.assign(nbperp, 0.0);
    if (Qsqr > 0) out.F_L_sqr.assign(nbperp, 0.0);
    out.F_T_integrand_sqr.assign(nbperp, 0.0);
    if (Qsqr > 0) out.F_L_integrand_sqr.assign(nbperp, 0.0);

    const double db = maxb / nbperp;
    for (int ib=0; ib<nbperp; ++ib)
        out.b[ib] = (ib + 0.5) * db;

    #pragma omp parallel for schedule(dynamic)
    for (int ib=0; ib<nbperp; ++ib) {
        const double bval = out.b[ib];
        // T polarization (vector integration returns real & imag)
        double int_modsq_T = 0.0;
        out.F_T[ib] = ScatteringAmplitudeF(xpom, Qsqr, bval, T, &int_modsq_T);
        out.F_T_sqr[ib] = std::norm(out.F_T[ib]);
        out.F_T_integrand_sqr[ib] = int_modsq_T;
        if (Qsqr > 0) {
            double int_modsq_L = 0.0;
            out.F_L[ib] = ScatteringAmplitudeF(xpom, Qsqr, bval, L, &int_modsq_L);
            out.F_L_sqr[ib] = std::norm(out.F_L[ib]);
            out.F_L_integrand_sqr[ib] = int_modsq_L;
        }
    }

    // Integrate over b for total cross sections in nb
    double sigma_T = 0.0, sigma_L = 0.0;
    for (int ib=0; ib<nbperp; ++ib) {
        double bval = out.b[ib];
        sigma_T += out.F_T_sqr[ib] * bval * db;
        if (Qsqr > 0) sigma_L += out.F_L_sqr[ib] * bval * db;
    }
    sigma_T *= 1e7 * HBARC * HBARC / (8. * M_PI);
    sigma_L *= 1e7 * HBARC * HBARC / (8. * M_PI);
    out.sigma_T = sigma_T;
    out.sigma_L = sigma_L;
    return out;
}

/*
 * Calculate scattering amplitude assuming that dipole amplitude does not depend on any angles
 * this is true if we have no constituent quarks in the ipsat model
 * No need to do mc integral, so this is numerically more accurate
 */
const int INTPOINTS_ROTSYM=2000;
double inthelperf_amplitude_rotationalsym_b(double b, void* p);
double inthelperf_amplitude_rotationalsym_r(double r, void* p);
double inthelperf_amplitude_rotationalsym_z(double z, void* p);

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
    double z;
    Polarization polarization;
};

double Diffraction::ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t, Polarization pol)
{
    
    Inthelper_amplitude par;
    par.diffraction = this;
    par.xpom=xpom; par.Qsqr = Qsqr; par.t=t;
    par.polarization= pol;
    
    gsl_function f;
    f.params = &par;
    f.function = inthelperf_amplitude_rotationalsym_b;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS_ROTSYM);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 100, 0, 0.001, INTPOINTS_ROTSYM, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#bint failed, result " << result << " relerror " << std::abs(error/result) << " t " <<t << endl;
    
    gsl_integration_workspace_free(w);
    
    if (std::isnan(result))
    {
        cerr<< "Diffraction::ScatteringAmplitudeRotationalSymmetry result is NaN, xpom=" << xpom << " t=" << t << endl;
    }
    
    return result;
}

double Diffraction::ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z, Polarization pol)
{
    // Set quark and antiquark on x axis around impact parameter
    // As amplitude does not depend on angle, this is ok here.
    Vec q1(b+r/2.0,0);
    Vec q2(b-r/2.0, 0);
    double amp = 2.0*dipole->Amplitude(xpom, q1, q2);
    //double overlap =wavef->PsiSqr_tot(Qsqr, r, z)/(4.0*M_PI);
    double overlap=0;
    if (pol == T){
        //overlap = wavef->PsiSqr_T_intz(Qsqr, r)/(4.0*M_PI);
        overlap = wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI);
    }
    else if (pol == L)
    {
        //overlap = wavef->PsiSqr_L_intz(Qsqr, r)/(4.0*M_PI);
        overlap = wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
    }
    else
        cerr << "Unknown polarization in Diffraction::ScatteringAmplitudeRotationalSymmetryIntegrand! " << endl;

    // Bessel integrals INCLUDING jacobian
    double delta = std::sqrt(t);
    double bessel = 2.0*M_PI*b*gsl_sf_bessel_J0(b*delta)*2.0*M_PI*r*gsl_sf_bessel_J0((1.0-z)*r*delta);
    
    return amp*overlap*bessel;
}

double inthelperf_amplitude_rotationalsym_b(double b, void* p)
{
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    par->b = b;
    gsl_function f;
    f.params = par;
    f.function = inthelperf_amplitude_rotationalsym_r;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS_ROTSYM);
    double result,error;
    int status = gsl_integration_qag(&f, std::log(1e-6), std::log(50), 0, 0.001, INTPOINTS_ROTSYM, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#R int failed, result " << result << " relerror " << std::abs(error/result) << " b " << b << " t " << par->t << endl;
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double inthelperf_amplitude_rotationalsym_r(double lnr, void* p)
{
    double r = exp(lnr);
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    par->r = r;
    // factorize z integral
    //return inthelperf_amplitude_rotationalsym_z(0.5, par);
    
    gsl_function f;
    f.params = par;
    f.function = inthelperf_amplitude_rotationalsym_z;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTPOINTS_ROTSYM);
    double result,error;
    double eps=1e-4;
    int status = gsl_integration_qags(&f, 0+eps, 1.0-eps, 0, 0.001, INTPOINTS_ROTSYM, w, &result, &error);
    
    
    if (status)
    {
        if (!(std::isnan(result)) and std::abs(result)>1e-20)
            cerr << "#zint failed, result " << result << " relerror " << std::abs(error/result) << " b " << par->b << " t " <<par->t << endl;
        if (std::isnan(result) and par->b < 30)
            cerr << " Nan also at b=" << par->b << endl;
        result=0;
    }
    
    
    gsl_integration_workspace_free(w);
    
    return r*result;    // r from exp(r) integration
    
}

double inthelperf_amplitude_rotationalsym_z(double z, void* p)
{
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    
    // Note: no jacobian here, it is included in ScatteringAmplitudeRotationalSymmetryIntegrand
    return par->diffraction->ScatteringAmplitudeRotationalSymmetryIntegrand(par->xpom, par->Qsqr, par->t, par->r, par->b, z, par->polarization);
}

// Helpers


DipoleAmplitude
* Diffraction::GetDipole()
{
    return dipole;
}
WaveFunction* Diffraction::GetWaveFunction(){
    return wavef;
}
