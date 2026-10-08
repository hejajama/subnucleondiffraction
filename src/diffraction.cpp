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

#include <algorithm>
#include <cmath>
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
        // Derived once per integral instead of at every integrand evaluation
        bool nrqcd;
        double umin;
        double umax;
        double delta;
        double J;
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
        factorize_zint,
        pol
    };
    p.nrqcd = wavef->WaveFunctionType() == "NRQCD";
    p.umin = std::log(p.rmin);
    p.umax = std::log(p.rmax);
    p.delta = std::sqrt(t);
    // Jacobian: bmax * (umax-umin) * (2pi)^2 * z-width (if unfactorized)
    const double twoPi = 2.0*M_PI;
    p.J = p.bmax * (p.umax-p.umin) * (twoPi*twoPi) * (p.fact ? 1.0 : (1.0 - 2.0*p.zmin));
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
        const double umin = prm->umin, umax = prm->umax;
        const double u = umin + (umax-umin) * xu;
        const double r = std::exp(u);
        const double theta_b = twoPi * xtb;
        const double theta_r = twoPi * xtr;
        const double z = fact ? 0.5 : (prm->zmin + (1.0 - 2.0*prm->zmin) * xz);
        // Overlap and scalar prefactor (2 r b * overlap)
        double scalar = 2.0 * r * b;
        if (fact) {
            if (prm->nrqcd) {
                if (prm->pol == T)
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_T_intz(prm->Q2, r, prm->delta, theta_r);
                else
                    scalar *= ((NRQCD_WF*)prm->diff->wavef)->PsiSqr_L_intz(prm->Q2, r, prm->delta, theta_r);
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
        
        
        // Quark coordinates
        // Note: as b is the center of the dipole, not the center-of-mass, no z factors here, 
        // but instead we have the off-forward phase below
        double qx, qy, qbarx, qbary;
        qx = bx + 0.5*rx; qy = by + 0.5*ry;
        qbarx = bx - 0.5*rx; qbary = by - 0.5*ry;

        double x1[2] = {qx,qy}; double x2[2] = {qbarx,qbary};
        std::complex<double> amp = prm->diff->dipole->ComplexAmplitude(prm->xp, x1, x2);
        // Phase factor with momentum transfer delta
        const double delta = prm->delta;
        if (delta > 0) {
            double phi = b*delta*cos_tb - (0.5 - z)*r*delta*cos_tr;

            // exp(-i phi)
            const std::complex<double> exponent(std::cos(phi), -std::sin(phi));
            amp *= exponent;
        }
        std::complex<double> val = scalar * amp;
        const double J = prm->J;
        const double measure_r = r; // from dr = r du
        f[0] = J * measure_r * static_cast<cubareal>(val.real());
        f[1] = J * measure_r * static_cast<cubareal>(val.imag());
        return 0;
    };
    const int ndim = factorize_zint ? 4 : 5;
    const int ncomp = 2;
    int nregions=0, neval=0, fail=0; double integral[2], error[2], prob[2];
    const int nvec = 1; const double epsrel = MCINTACCURACY, epsabs = 0.0;
    const int flags = 0, seed = 0;
    const int mineval = mcintpoints/10; const int maxeval = mcintpoints;
    // Suave fails to allocate its regions if nnew < nmin (mcintpoints < 60000);
    // it then still evaluates at least nmin points (main() rejects fewer)
    const int nmin = 300; const int nnew = std::max(mineval/20, nmin); const double flatness = 1.0;
    Suave(ndim, ncomp, integrand, &p, nvec, epsrel, epsabs, flags, seed,
        mineval, maxeval, nnew, nmin, flatness,
        NULL, NULL, &nregions, &neval, &fail, integral, error, prob);
    return std::complex<double>(integral[0], integral[1]);
}

std::complex<double> Diffraction::ScatteringAmplitude_tIntegrated(
    double xpom, double Qsqr, double b, double theta_b, Polarization pol, double epsabs) {
    struct SuaveParams {
        Diffraction* diff;
        double xpom;
        double Q2;
        double b;
        double theta_b;
        double zmin;
        double rmin;
        double rmax;
        bool factorize;
        Polarization pol;
        // Derived once per integral instead of at every integrand evaluation
        bool nrqcd;
        double umin;
        double umax;
        double J;
        long nonfinite; // number of integrand evaluations that returned NaN/inf
    };
    SuaveParams p{
        this,
        xpom,
        Qsqr,
        b,
        theta_b,
        zlimit,
        1e-10,
        MAXR,
        factorize_zint,
        pol
    };
    p.nrqcd = wavef->WaveFunctionType() == "NRQCD";
    p.umin = std::log(p.rmin);
    p.umax = std::log(p.rmax);
    // Overall Jacobian (theta_r, u, z (if not factorized))
    p.J = 2.0*M_PI * (p.umax-p.umin) * (p.factorize ? 1.0 : (1.0 - 2.0*p.zmin));
    p.nonfinite = 0;
    auto integrand = [](const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *ud)->int{
        SuaveParams* prm = static_cast<SuaveParams*>(ud);
        const double twoPi = 2.0*M_PI;
        const bool fact = prm->factorize;


        const double umin = prm->umin, umax = prm->umax;
        const double xr = x[0];
        const double xu = x[1];
        const double xz = fact ? 0.5 : x[2];
        const double theta_b_int = prm->theta_b;
        const double theta_r = twoPi * xr;
        const double u = umin + (umax-umin) * xu;
        const double r = std::exp(u);
        const double z = fact ? 0.5 : (prm->zmin + (1.0 - 2.0*prm->zmin) * xz);
        // Common factors
        double scalar = r; // r from Jacobian (du->dr adds r)
        if (fact) {
            if (prm->nrqcd) {
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
            if (prm->pol == T)
                scalar *= prm->diff->wavef->PsiSqr_T(prm->Q2, r, z);
            else
                scalar *= prm->diff->wavef->PsiSqr_L(prm->Q2, r, z);
        }
       
        const double cos_tb = std::cos(theta_b_int);
        const double sin_tb = std::sin(theta_b_int);
        const double cos_tr = std::cos(theta_r);
        const double sin_tr = std::sin(theta_r);
        double bx = prm->b * cos_tb;
        double by = prm->b * sin_tb;
        double rx = r * cos_tr;
        double ry = r * sin_tr;
        // q and antiq positions
        // Note: no off forward phase, instead b is the center-of-mass of the dipole
        double qx = bx + (1. - z) * rx;
        double qy = by + (1. - z) * ry;
        double qbarx = bx - z * rx;
        double qbary = by - z * ry;
        double x1[2] = {qx,qy};
        double x2[2] = {qbarx,qbary};
        std::complex<double> amp = prm->diff->dipole->ComplexAmplitude(prm->xpom, x1, x2);
        const double amp_r = amp.real();
        const double amp_i = amp.imag();
        if (!std::isfinite(amp_r) || !std::isfinite(amp_i) || !std::isfinite(scalar)) {
            // Do not let a single bad point poison the Suave variance estimate
            ++prm->nonfinite;
            f[0] = 0.0; f[1] = 0.0;
            return 0;
        }
        const double J = prm->J;
        // Jacobian pieces:
        //  theta_r: 2pi  (in J)
        //  u = ln r mapping: u = umin + (umax-umin)*xu gives width (umax-umin) in J and dr = r du adds extra r
        //  optional z: width (1 - 2 zmin) in J
        // scalar currently includes factor r * overlap; needs extra r from dr=r du
        const double measure_u_r = r; // u = ln r
        // Components: real, imag
        f[0] = J * measure_u_r * scalar * amp_r; // real part
        f[1] = J * measure_u_r * scalar * amp_i; // imag part
        return 0;
    };

    if (factorize_zint)
    {   
        cerr << "factorize_zint in ScatteringAmpltiude_tIntegrated has not been tested" << endl;
        exit(1);
    }

    const int ndim = factorize_zint ? 2 : 3;
    const int ncomp = 2; // real, imag
    int nregions=0, neval=0, fail=0;
    double integral[2], error[2], prob[2];
    const int nvec = 1;
    const double epsrel = MCINTACCURACY;
    const int flags = 0, seed = 0;
    const int mineval = mcintpoints/10; const int maxeval = mcintpoints;
    // Suave fails to allocate its regions if nnew < nmin (mcintpoints < 60000);
    // it then still evaluates at least nmin points (main() rejects fewer)
    const int nmin = 300; const int nnew = std::max(mineval/20, nmin); const double flatness = 1.0;

    Suave(ndim, ncomp, integrand, &p, nvec, epsrel, epsabs, flags, seed,
        mineval, maxeval, nnew, nmin, flatness,
        NULL, NULL, &nregions, &neval, &fail, integral, error, prob);
    // fail != 0 alone only means that the relative accuracy goal was not
    // reached within maxeval, which is the normal case here
    if (p.nonfinite > 0 || !std::isfinite(integral[0]) || !std::isfinite(integral[1])) {
        #pragma omp critical
        cerr << "# ScatteringAmplitude_tIntegrated: Suave fail=" << fail << " b=" << b
             << " theta_b=" << theta_b << " neval=" << neval << " nregions=" << nregions
             << " nonfinite_evals=" << p.nonfinite
             << " result=(" << integral[0] << " +- " << error[0] << ", "
             << integral[1] << " +- " << error[1] << ")" << endl;
    }
    return std::complex<double>(integral[0], integral[1]);
}

Diffraction::TotalCrossSectionData Diffraction::ComputeTotalCrossSection(
    double xpom, double Qsqr, int nbperp, double maxb, int ntheta) {
    TotalCrossSectionData out;
    out.b.resize(nbperp);
    out.theta.resize(ntheta);
    const int ntot = nbperp * ntheta;
    out.F_T.assign(ntot, std::complex<double>(0.,0.));
    if (Qsqr > 0) out.F_L.assign(ntot, std::complex<double>(0.,0.));

    const double db = maxb / nbperp;
    for (int ib=0; ib<nbperp; ++ib)
        out.b[ib] = (ib + 0.5) * db;

    const double dtheta = 2.0 * M_PI / ntheta;
    for (int it=0; it<ntheta; ++it)
        out.theta[it] = it * dtheta;

    // With epsabs = 0, Suave uses all maxeval points also where the integrand is
    // (practically) zero, at large b. Set an absolute accuracy goal relative to
    // the amplitude at the smallest b, so that it stops early there; elsewhere
    // the error stays far above it and the results are unchanged.
    std::complex<double> F_ref = ScatteringAmplitude_tIntegrated(xpom, Qsqr, out.b[0], 0.0, T);
    double scale = std::abs(F_ref);
    if (Qsqr > 0)
        scale = std::max(scale, std::abs(ScatteringAmplitude_tIntegrated(xpom, Qsqr, out.b[0], 0.0, L)));
    const double epsabs = 1e-3 * MCINTACCURACY * scale;

    #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int ib=0; ib<nbperp; ++ib) {
        for (int it=0; it<ntheta; ++it) {
            const int idx = ib*ntheta + it;
            const double bval = out.b[ib];
            const double thetaval = out.theta[it];
            // T polarization (vector integration returns real & imag)
            out.F_T[idx] = ScatteringAmplitude_tIntegrated(xpom, Qsqr, bval, thetaval, T, epsabs);
            if (Qsqr > 0) {
                out.F_L[idx] = ScatteringAmplitude_tIntegrated(xpom, Qsqr, bval, thetaval, L, epsabs);
            }
        }
    }

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
