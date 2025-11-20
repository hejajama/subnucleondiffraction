/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015-2016
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
#include <omp.h>
#include <functional>
#include <iomanip>


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    num_of_averages = 1;
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
    bool real_part;
    Polarization polarization;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

double Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, Polarization pol, bool real_part)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization = pol;
    helper.real_part = real_part;

    // Do MC integral over impact parameters and dipole sizes    
    // Currently hardcoded parameters for jpsi and gold:
    // Impact parameter up to 100 GeV^-1
    // Dipole size up to 10 GeV^-1
    // MC integration parameters: b, theta_b, r, theta_r, Z
    double *lower, *upper;
    if (FACTORIZE_ZINT)
    {
        lower = new double[4];
        upper = new double[4];
    }
    else
    {
        lower = new double[5];
        upper = new double[5];
        lower[4] = zlimit; // Min z
        upper[4] = 1.0 - lower[4];    // Max z
    }

    lower[0] = lower[1] = lower[2] = lower[3] = 0;
    upper[0] = 10 * 5.068; // Max b
    upper[1] = 2.0 * M_PI;
    upper[2] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[3] = 2.0 * M_PI;

    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 4;
    if (!FACTORIZE_ZINT)
        F.dim = 5;
    F.params = &helper;

    double result = 0.;
    double error = 0.;

    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS,
                                  global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    } else if (MCINT == VEGAS) {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50,
                                  global_rng, s, &result, &error);

        if (ShowVegasIterations() == false)
            cerr << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5,
                                      global_rng, s, &result, &error);
            if (ShowVegasIterations() == false)
                cerr << "# Vegas integration " << result << " +/- " << error
                     << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            if (iter>10)
                break;
        } while (iter < 2 or ( std::abs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 or std::abs(error/result) > MCINTACCURACY));
        gsl_monte_vegas_free(s);
    }
    delete[] lower;
    delete[] upper;
    return result;
}


double Inthelperf_amplitudeF_mc(double *vec, size_t dim, void* par);

// Deterministic (GSL) helpers for ScatteringAmplitudeF
struct FDetParams {
    Diffraction* diff;
    double xpom;
    double Q2;
    double b;
    double theta_b;
    double theta_r;
    Polarization pol;
    bool real_part;
    double zmin;
    gsl_integration_workspace* wz;
};

struct ZIntParams { FDetParams* base; double r; };

static double FDetZHelper(double z, void* p)
{
    ZIntParams* prm = static_cast<ZIntParams*>(p);
    FDetParams* b = prm->base;
    return b->diff->ScatteringAmplitudeIntegrand2(
        b->xpom, b->Q2, prm->r, b->theta_r, b->b, b->theta_b, z, b->pol, b->real_part);
}

static double FDetRLogHelper(double u, void* p)
{
    FDetParams* prm = static_cast<FDetParams*>(p);
    const double r = std::exp(u);
    if (FACTORIZE_ZINT) {
        double f = prm->diff->ScatteringAmplitudeIntegrand2(
            prm->xpom, prm->Q2, r, prm->theta_r, prm->b, prm->theta_b, 0.5, prm->pol, prm->real_part);
        return r * f; // dr/du = r
    } else {
        ZIntParams zp; zp.base = prm; zp.r = r;
        double zl = std::max(0.0, prm->zmin);
        double zh = std::min(1.0, 1.0 - prm->zmin);

        // First try fast fixed Gauss–Legendre with 16 points (and 8 for error estimate)
        static thread_local gsl_integration_glfixed_table* gl_z_16 = nullptr;
        static thread_local gsl_integration_glfixed_table* gl_z_8 = nullptr;
        if (!gl_z_16) gl_z_16 = gsl_integration_glfixed_table_alloc(16);
        if (!gl_z_8) gl_z_8 = gsl_integration_glfixed_table_alloc(8);
        auto gl_eval = [&](gsl_integration_glfixed_table* tab)->double {
            double acc = 0.0;
            for (int i=0; i<tab->n; ++i) {
                double zi, wi;
                gsl_integration_glfixed_point(zl, zh, i, &zi, &wi, tab);
                acc += wi * FDetZHelper(zi, &zp);
            }
            return acc;
        };
        double z16 = gl_eval(gl_z_16);
        double z8  = gl_eval(gl_z_8);
        double denom = std::max(1.0, std::abs(z16));
        if (std::abs(z16 - z8) / denom <= ZINT_RELACCURACY) {
            return r * z16;
        }

        // If GL not sufficient, try qng then fallback to qag
        gsl_function Fz; Fz.function = FDetZHelper; Fz.params = &zp;
        double zres=0.0, zerr=0.0; size_t neval = 0;
        int st = gsl_integration_qng(&Fz, zl, zh, 0.0, ZINT_RELACCURACY, &zres, &zerr, &neval);
        if (st != GSL_SUCCESS) {
            st = gsl_integration_qag(&Fz, zl, zh, 0.0, ZINT_RELACCURACY, ZINT_INTERVALS, GSL_INTEG_GAUSS51, prm->wz, &zres, &zerr);
            if (st) {
                if (std::isnan(zres) || std::isinf(zres)) return 0.0;
            }
        }
        return r * zres;
    }
}

double Diffraction::ScatteringAmplitudeF(
    double xpom, double Qsqr, double b, double theta_b, Polarization pol, bool real_part) {

    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = 0.0;
    helper.b = b;
    helper.theta_b = theta_b;
    helper.polarization = pol;
    helper.real_part = real_part;

    // Do MC integral over impact parameters and dipole sizes
    // Currently hardcoded parameters for jpsi and gold:
    // Impact parameter up to 100 GeV^-1
    // Dipole size up to 10 GeV^-1
    // MC integration parameters: r, theta_r, z
    double *lower, *upper;
    if (FACTORIZE_ZINT)
    {
        lower = new double[2];
        upper = new double[2];
    }
    else
    {
        lower = new double[3];
        upper = new double[3];
        lower[2] = zlimit; // Min z
        upper[2] = 1. - lower[2];    // Max z
    }
    lower[0] = lower[1] = 0;
    upper[0] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[1] = 2.0*M_PI;        // phi_r

    gsl_monte_function F;
    F.f = &Inthelperf_amplitudeF_mc;
    F.dim = 2;
    if (!FACTORIZE_ZINT)
        F.dim = 3;
    F.params = &helper;

    double result = 0.;
    double error = 0.;

    if (MCINT == MISER) {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS,
                                  global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    } else if (MCINT == VEGAS) {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50,
                                  global_rng, s, &result, &error);

        if (ShowVegasIterations())
            cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5,
                                      global_rng, s, &result, &error);
            if (ShowVegasIterations())
                cout << "# Vegas integration " << result << " +/- " << error
                     << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            if (iter>10)
                break;
        } while (iter < 2 or ( std::abs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 or std::abs(error/result) > MCINTACCURACY));
        gsl_monte_vegas_free(s);
    } else if (MCINT == GSL) {
        // Deterministic integration with adaptive composite Simpson in theta_r:
        // 1. Start with 32 base intervals (uniform partition).
        // 2. For each interval compute quarter-point refinement only once; cache error and refined Simpson estimate.
        // 3. Apply Richardson correction ( (S_refined - S_base)/15 ) when accepting interval.
        // 4. Refine (split) only intervals whose relative error > MCINTACCURACY, until cap of theta evaluations.
        static thread_local gsl_integration_cquad_workspace* wrcq = gsl_integration_cquad_workspace_alloc(1024);

        FDetParams prm; prm.diff = this; prm.xpom = xpom; prm.Q2 = Qsqr; prm.b = b; prm.theta_b = theta_b;
        prm.pol = pol; prm.real_part = real_part; prm.zmin = zlimit; prm.wz = gsl_integration_workspace_alloc(ZINT_INTERVALS);

        auto r_integrated = [&](double theta_r)->double {
            prm.theta_r = theta_r;
            gsl_function Fr; Fr.function = FDetRLogHelper; Fr.params = &prm;
            const double rmin = 1e-6; const double a = std::log(rmin); const double bb = std::log(MaxR());
            double rres = 0.0, rerr = 0.0; size_t neval = 0;
            int st = gsl_integration_qng(&Fr, a, bb, 0.0, MCINTACCURACY, &rres, &rerr, &neval);
            if (st != GSL_SUCCESS || !(std::abs(rres) > 0 && std::abs(rerr/rres) <= MCINTACCURACY)) {
                st = gsl_integration_cquad(&Fr, a, bb, 0.0, MCINTACCURACY, wrcq, &rres, &rerr, &neval);
                if (st && (std::isnan(rres) || std::isinf(rres))) rres = 0.0;
            }
            return rres;
        };

        struct ThetaInterval {
            double a, b; // bounds
            double fa, fm, fb; // values at a, midpoint, b
            double Sbase;      // base Simpson over [a,b]
            // Cached refinement data
            bool error_computed;
            double fql, fqr;   // quarter points
            double Srefined;   // refined Simpson (sum of two halves)
            double relerr;     // relative error estimate
        };

        const double A = 0.0, B = 2.0*M_PI;
        const int baseIntervals = 16;
        const int maxThetaEvalsCap = 256;
        std::vector<ThetaInterval> intervals; intervals.reserve(baseIntervals*2);
        const double h = (B - A) / baseIntervals;

        int evals = 0;
        // Precompute f(0) exploiting periodicity later
        double f0 = r_integrated(A); evals++;
        for (int i=0; i<baseIntervals; ++i) {
            double a_i = A + i*h;
            double b_i = a_i + h;
            double fa = (i==0 ? f0 : r_integrated(a_i)); if (i>0) evals++;
            double fb = (i==baseIntervals-1 ? f0 : r_integrated(b_i)); if (i!=baseIntervals-1) evals++;
            double m  = 0.5*(a_i + b_i);
            double fm = r_integrated(m); evals++;
            double Sbase = (h/6.0)*(fa + 4.0*fm + fb);
            intervals.push_back({a_i,b_i,fa,fm,fb,Sbase,false,0.0,0.0,0.0,0.0});
        }

        auto compute_error = [&](ThetaInterval& I){
            if (I.error_computed) return; // already done
            if (evals >= maxThetaEvalsCap) return; // cannot refine further
            double m = 0.5*(I.a + I.b);
            double leftMid  = 0.5*(I.a + m);
            double rightMid = 0.5*(m + I.b);
            I.fql = r_integrated(leftMid); evals++;
            if (evals < maxThetaEvalsCap) { I.fqr = r_integrated(rightMid); evals++; }
            else I.fqr = I.fql; // fallback if cap reached (degenerate refinement)
            double Sl = ( (m - I.a)/6.0 ) * ( I.fa + 4.0*I.fql + I.fm );
            double Sr = ( (I.b - m)/6.0 ) * ( I.fm + 4.0*I.fqr + I.fb );
            I.Srefined = Sl + Sr;
            double denom = std::max(1.0, std::abs(I.Srefined));
            double absErr = std::abs(I.Srefined - I.Sbase) / 15.0; // classical Simpson error numerator
            I.relerr = absErr / denom;
            I.error_computed = true;
        };

        // Initial one-pass error computation for all intervals; then assess GLOBAL errors
        for (auto& I : intervals) compute_error(I);
        double sum_base = 0.0, sum_refined = 0.0, sum_abs_diff = 0.0;
        double worstRelErr = 0.0;
        for (auto& I : intervals) {
            sum_base += I.Sbase;
            sum_refined += I.Srefined;
            sum_abs_diff += std::abs(I.Srefined - I.Sbase);
            if (I.relerr > worstRelErr) worstRelErr = I.relerr;
        }
        double denom_global = std::max(1.0, std::abs(sum_refined));
        double global_abs_err = (sum_abs_diff / denom_global) / 15.0; // absolute-sum composite Simpson error

        int splitsPerformed = 0;
        // Hybrid refinement: continue while either worst per-interval or global abs error exceed tolerance
        while (evals < maxThetaEvalsCap && (worstRelErr > MCINTACCURACY || global_abs_err > MCINTACCURACY)) {
            // Pick worst interval by relerr
            double worstErrLocal = -1.0; size_t worstIdx = intervals.size();
            for (size_t i=0; i<intervals.size(); ++i) {
                if (intervals[i].relerr > worstErrLocal) { worstErrLocal = intervals[i].relerr; worstIdx = i; }
            }
            if (worstIdx == intervals.size()) break; // no candidate
            ThetaInterval Iold = intervals[worstIdx];
            if (evals >= maxThetaEvalsCap) break;
            double m = 0.5*(Iold.a + Iold.b);
            // Children (use cached quarter points as midpoints)
            ThetaInterval left{Iold.a, m, Iold.fa, Iold.fql, Iold.fm,
                               ( (m - Iold.a)/6.0 ) * ( Iold.fa + 4.0*Iold.fql + Iold.fm ), false,0,0,0,0};
            ThetaInterval right{m, Iold.b, Iold.fm, Iold.fqr, Iold.fb,
                                ( (Iold.b - m)/6.0 ) * ( Iold.fm + 4.0*Iold.fqr + Iold.fb ), false,0,0,0,0};
            intervals[worstIdx] = left;
            intervals.insert(intervals.begin() + worstIdx + 1, right);
            // Compute errors for new children
            compute_error(intervals[worstIdx]);
            compute_error(intervals[worstIdx+1]);
            // Recompute global metrics
            sum_base = 0.0; sum_refined = 0.0; sum_abs_diff = 0.0; worstRelErr = 0.0;
            for (auto& I : intervals) {
                sum_base += I.Sbase; sum_refined += I.Srefined;
                sum_abs_diff += std::abs(I.Srefined - I.Sbase);
                if (I.relerr > worstRelErr) worstRelErr = I.relerr;
            }
            denom_global = std::max(1.0, std::abs(sum_refined));
            global_abs_err = (sum_abs_diff / denom_global) / 15.0;
            splitsPerformed++;
        }

        // Final accumulation with Richardson correction on refined intervals
        double sum = 0.0; int acceptedRefined = 0; int totalIntervals = (int)intervals.size();
        for (auto& I : intervals) {
            double contrib = I.Sbase; // default base Simpson
            if (I.error_computed) { contrib = I.Srefined + (I.Srefined - I.Sbase)/15.0; acceptedRefined++; }
            sum += contrib;
        }
        result = sum;
        /*cerr << "# GSL theta_r evals " << evals << " intervals " << totalIntervals
             << " refined_passes " << acceptedRefined
             << " splits " << splitsPerformed
             << " worst_relerr " << worstRelErr
             << " global_abs_err " << global_abs_err
             << " (cap " << maxThetaEvalsCap << ") for b=" << b << endl;*/
        gsl_integration_workspace_free(prm.wz);
    }
    delete[] lower;
    delete[] upper;
    return result;
}

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->b = vec[0];
    helper->theta_b=vec[1];
    helper->r = vec[2];
    helper->theta_r = vec[3];
    
    double z = 0.5;// Put z=0.5 as it sets b to the geometric average of quarks
    
    if (!FACTORIZE_ZINT)
        z = vec[4];
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(
        helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, 
        helper->b, helper->theta_b, z, helper->polarization, helper->real_part);
}

double Inthelperf_amplitudeF_mc(double *vec, size_t dim, void* par) {
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->r = vec[0];
    helper->theta_r = vec[1];

    double z = 0.5;// Put z=0.5 as it sets b to the geometric average of quarks

    if (!FACTORIZE_ZINT)
        z = vec[2];

    return helper->diffraction->ScatteringAmplitudeIntegrand2(
        helper->xpom, helper->Qsqr, helper->r, helper->theta_r,
        helper->b, helper->theta_b, z, helper->polarization, helper->real_part);
}

double Diffraction::ScatteringAmplitudeIntegrand(
    double xpom, double Qsqr, double t, double r, double theta_r,
    double b, double theta_b, double z, Polarization pol, bool real_part) {

    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    // If do like Lappi, Mantysaari: set z=1/2 here

    const double cos_theta_b = std::cos(theta_b);
    const double sin_theta_b = std::sin(theta_b);
    const double cos_theta_r = std::cos(theta_r);
    const double sin_theta_r = std::sin(theta_r);

    double bx = b*cos_theta_b;
    double by = b*sin_theta_b;
    double rx = r*cos_theta_r;
    double ry = r*sin_theta_r;

    // q and antiq positions
    // Note my convention is that b is the center of the dipole (geometric center), not center of mass (z weighted)
    // Consequently I get (0.5-z)r.Delta phase    
    double qx = bx + 0.5*rx; double qy = by + 0.5*ry;
    double qbarx = bx - 0.5*rx; double qbary = by - 0.5*ry;

    double delta = std::sqrt(t);

    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    std::complex<double> amp = dipole->ComplexAmplitude(xpom, x1, x2);
    const double amp_r = amp.real();
    const double amp_i = amp.imag();
    double scalar = 2.0*r*b; // r and b from Jacobians, 2 as we have written sigma_qq = 2 N

    double overlap = 0.0;
    if (FACTORIZE_ZINT) {
        if (wavef->WaveFunctionType() == "NRQCD") {
            if (pol == T)
                overlap = ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
            else
                overlap = ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta, theta_r);
        } else {
            if (pol == T)
                overlap = wavef->PsiSqr_T_intz(Qsqr, r);
            else
                overlap = wavef->PsiSqr_L_intz(Qsqr, r);
        }
        scalar *= overlap;
        double res;
        if (delta > 0) {
            const double phi = - b*delta*cos_theta_b;
            const double c = std::cos(phi);
            const double s = std::sin(phi);
            if (real_part)
                res = scalar * (amp_r*c - amp_i*s);
            else
                res = scalar * (amp_r*s + amp_i*c);
        } else {
            res = real_part ? scalar * amp_r : scalar * amp_i;
        }
        if (std::isnan(res) || std::isinf(res)) {
            cerr << "Amplitude integral is " << res << " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
            exit(1);
        }
        return res;
    } else {
        const double inv4pi = 1.0/(4.0*M_PI);
        if (pol == T)
            overlap = wavef->PsiSqr_T(Qsqr, r, z) * inv4pi;
        else
            overlap = wavef->PsiSqr_L(Qsqr, r, z) * inv4pi;
        scalar *= overlap;
        double res;
        if (delta > 0) {
            const double phi = - ( b*delta*cos_theta_b - (0.5 - z)*r*delta*cos_theta_r );
            const double c = std::cos(phi);
            const double s = std::sin(phi);
            if (real_part)
                res = scalar * (amp_r*c - amp_i*s);
            else
                res = scalar * (amp_r*s + amp_i*c);
        } else {
            res = real_part ? scalar * amp_r : scalar * amp_i;
        }
        if (std::isnan(res) || std::isinf(res)) {
            cerr << "Amplitude integral is " << res << " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
            exit(1);
        }
        return res;
    }
}


double Diffraction::ScatteringAmplitudeIntegrand2(
    double xpom, double Qsqr, double r, double theta_r,
    double b, double theta_b, double z, Polarization pol, bool real_part) { 

    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    // If do like Lappi, Mantysaari: set z=1/2 here

    const double cos_theta_b = std::cos(theta_b);
    const double sin_theta_b = std::sin(theta_b);
    const double cos_theta_r = std::cos(theta_r);
    const double sin_theta_r = std::sin(theta_r);

    double bx = b*cos_theta_b;
    double by = b*sin_theta_b;
    double rx = r*cos_theta_r;
    double ry = r*sin_theta_r;

    // q and antiq positions
    double qx = bx + (1. - z) * rx;
    double qy = by + (1. - z) * ry;
    double qbarx = bx - z * rx;
    double qbary = by - z * ry;

    double x1[2] = {qx, qy};
    double x2[2] = {qbarx, qbary};
    std::complex<double> amp = dipole->ComplexAmplitude(xpom, x1, x2);
    const double amp_r = amp.real();
    const double amp_i = amp.imag();
    double scalar = 2.0*r; // r from Jacobians, 2 as we have written sigma_qq = 2 N

    if (FACTORIZE_ZINT) {
        if (wavef->WaveFunctionType() == "NRQCD") {
            double delta = 0.0;
            if (pol == T)
                scalar *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
            else
                scalar *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta, theta_r);
        } else {
            if (pol == T)
                scalar *= wavef->PsiSqr_T_intz(Qsqr, r);
            else
                scalar *= wavef->PsiSqr_L_intz(Qsqr, r);
        }
    } else {
        const double inv4pi = 1.0/(4.0*M_PI);
        if (pol == T)
            scalar *= wavef->PsiSqr_T(Qsqr, r, z) * inv4pi;
        else
            scalar *= wavef->PsiSqr_L(Qsqr, r, z) * inv4pi;
    }

    return real_part ? scalar * amp_r : scalar * amp_i;
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
    //double bessel=2.0*M_PI*b*2.0*M_PI*r*1.0;
    
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
    
    //cout << "rint at b=" << b << ": " << result << " pm " << error << endl;
    
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
    
    //cout << "zint at r=" << r << ", b=" << par->b << ": " << result << " pm " << error << endl;
    
    return r*result;    // r from exp(r) integration
    
}

double inthelperf_amplitude_rotationalsym_z(double z, void* p)
{
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    
    // Note: no jacobian here, it is included in ScatteringAmplitudeRotationalSymmetryIntegrand
    return par->diffraction->ScatteringAmplitudeRotationalSymmetryIntegrand(par->xpom, par->Qsqr, par->t, par->r, par->b, z, par->polarization);
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

Diffraction::TotalCrossSectionData Diffraction::ComputeTotalCrossSection(
    double xpom, double Qsqr, int nbperp, double maxb) {
    TotalCrossSectionData out;
    out.b.resize(nbperp);
    out.F_T.assign(nbperp, std::complex<double>(0.,0.));
    if (Qsqr > 0) out.F_L.assign(nbperp, std::complex<double>(0.,0.));
    out.F_T_sqr.assign(nbperp, 0.0);
    if (Qsqr > 0) out.F_L_sqr.assign(nbperp, 0.0);

    const double db = maxb / nbperp;
    for (int ib=0; ib<nbperp; ++ib)
        out.b[ib] = (ib + 0.5) * db;

    // Locally adaptive theta_b integration: start from a uniform base grid before refinement.
    const int maxEvalsCap = 256;          // hard cap on theta evaluations per b
    const int baseIntervalsThetaB = 16;   // minimum starting intervals (>= 3)

    #pragma omp parallel for schedule(dynamic)
    for (int ib=0; ib<nbperp; ++ib) {
        const double bval = out.b[ib];
        #pragma omp critical
        {
            cerr << "# Computing total cross section at b=" << bval << endl;
        }

        struct ThetaEval { std::complex<double> aT; std::complex<double> aL; double t2; double l2; };
        auto eval = [&](double th)->ThetaEval {
            double rT = ScatteringAmplitudeF(xpom, Qsqr, bval, th, T, true);
            double iT = ScatteringAmplitudeF(xpom, Qsqr, bval, th, T, false);
            std::complex<double> aT(rT, iT);
            double rL=0.0,iL=0.0; std::complex<double> aL(0.,0.);
            if (Qsqr > 0) {
                rL = ScatteringAmplitudeF(xpom, Qsqr, bval, th, L, true);
                iL = ScatteringAmplitudeF(xpom, Qsqr, bval, th, L, false);
                aL = std::complex<double>(rL, iL);
            }
            return ThetaEval{aT, aL, std::norm(aT), (Qsqr>0? std::norm(aL): 0.0)};
        };

        struct Interval {
            double a, b;          // bounds
            ThetaEval fa, fm, fb; // samples at a, mid, b
            // Base Simpson contributions
            std::complex<double> S_T; double S_T2; std::complex<double> S_L; double S_L2;
            // Refinement cache
            bool refined; ThetaEval fql, fqr; // quarter points
            std::complex<double> Sref_T; double Sref_T2; std::complex<double> Sref_L; double Sref_L2; double relerr;
        };
        struct SimpsonRes { std::complex<double> ST; double ST2; std::complex<double> SL; double SL2; };
        auto simpson = [&](const ThetaEval& fa, const ThetaEval& fm, const ThetaEval& fb, double h)->SimpsonRes {
            SimpsonRes R; R.ST = (h/6.0) * (fa.aT + 4.0*fm.aT + fb.aT);
            R.ST2 = (h/6.0) * (fa.t2 + 4.0*fm.t2 + fb.t2);
            R.SL = std::complex<double>(0.,0.); R.SL2 = 0.0;
            if (Qsqr > 0) {
                R.SL  = (h/6.0) * (fa.aL + 4.0*fm.aL + fb.aL);
                R.SL2 = (h/6.0) * (fa.l2 + 4.0*fm.l2 + fb.l2);
            }
            return R;
        };

        std::vector<Interval> intervals; intervals.reserve(baseIntervalsThetaB*2);
        const double A = 0.0, B = 2.0*M_PI;
        const double h0 = (B - A) / baseIntervalsThetaB;
        int evals = 0;
        // periodic value at 0 and 2pi
        ThetaEval f0 = eval(A); evals++;
        for (int i=0; i<baseIntervalsThetaB; ++i) {
            double a_i = A + i*h0;
            double b_i = a_i + h0;
            ThetaEval fa = (i==0? f0 : eval(a_i)); if (i>0) evals++;
            ThetaEval fb = (i==baseIntervalsThetaB-1? f0 : eval(b_i)); if (i!=baseIntervalsThetaB-1) evals++;
            double m = 0.5*(a_i + b_i);
            ThetaEval fm = eval(m); evals++;
            SimpsonRes SR = simpson(fa,fm,fb,h0);
            intervals.push_back(Interval{a_i,b_i,fa,fm,fb,SR.ST,SR.ST2,SR.SL,SR.SL2,false,ThetaEval(),ThetaEval(),std::complex<double>(0.,0.),0.0,std::complex<double>(0.,0.),0.0,0.0});
        }

        auto compute_error = [&](Interval& I){
            if (I.refined || evals >= maxEvalsCap) return;
            double m = 0.5*(I.a + I.b);
            double leftMid  = 0.5*(I.a + m);
            double rightMid = 0.5*(m + I.b);
            I.fql = eval(leftMid); evals++;
            if (evals < maxEvalsCap) { I.fqr = eval(rightMid); evals++; }
            else I.fqr = I.fql; // degenerate if cap reached
            // left Simpson
            SimpsonRes SRl = simpson(I.fa, I.fql, I.fm, m - I.a);
            SimpsonRes SRr = simpson(I.fm, I.fqr, I.fb, I.b - m);
            I.Sref_T  = SRl.ST + SRr.ST;
            I.Sref_T2 = SRl.ST2 + SRr.ST2;
            I.Sref_L  = SRl.SL + SRr.SL;
            I.Sref_L2 = SRl.SL2 + SRr.SL2;
            double denom = std::max(1.0, std::abs(I.Sref_T));
            double absErr = std::abs(I.Sref_T - I.S_T) / 15.0;
            I.relerr = absErr / denom;
            I.refined = true;
        };

        // Compute initial errors for all base intervals
        for (auto& I : intervals) if (!I.refined) compute_error(I);

        // Assess GLOBAL errors (absolute-sum to avoid cancellation) and worst interval error
        std::complex<double> sum_base_T(0.,0.), sum_refined_T(0.,0.);
        double sum_abs_diff_T = 0.0; double worst_relerr_T = 0.0;
        for (auto& I : intervals) {
            sum_base_T += I.S_T; sum_refined_T += I.Sref_T;
            sum_abs_diff_T += std::abs(I.Sref_T - I.S_T);
            if (I.relerr > worst_relerr_T) worst_relerr_T = I.relerr;
        }
        double denom_global_T = std::max(1.0, std::abs(sum_refined_T));
        double global_abs_err_T = (sum_abs_diff_T / denom_global_T) / 15.0;

        int splitsPerformedThetaB = 0;
        // Hybrid refinement: split while either worst per-interval or global absolute error exceed tolerance
        while (evals < maxEvalsCap && (worst_relerr_T > MCINTACCURACY || global_abs_err_T > MCINTACCURACY)) {
            double worstErrLocal = -1.0; size_t worstIdx = intervals.size();
            for (size_t i=0; i<intervals.size(); ++i) {
                if (intervals[i].relerr > worstErrLocal) { worstErrLocal = intervals[i].relerr; worstIdx = i; }
            }
            if (worstIdx == intervals.size()) break; // nothing to split
            if (evals >= maxEvalsCap) break;

            Interval Iold = intervals[worstIdx];
            double m = 0.5*(Iold.a + Iold.b);
            SimpsonRes SRl = simpson(Iold.fa, Iold.fql, Iold.fm, m - Iold.a);
            Interval left{Iold.a, m, Iold.fa, Iold.fql, Iold.fm, SRl.ST, SRl.ST2, SRl.SL, SRl.SL2, false,ThetaEval(),ThetaEval(),std::complex<double>(0.,0.),0.0,std::complex<double>(0.,0.),0.0,0.0};
            SimpsonRes SRr = simpson(Iold.fm, Iold.fqr, Iold.fb, Iold.b - m);
            Interval right{m, Iold.b, Iold.fm, Iold.fqr, Iold.fb, SRr.ST, SRr.ST2, SRr.SL, SRr.SL2, false,ThetaEval(),ThetaEval(),std::complex<double>(0.,0.),0.0,std::complex<double>(0.,0.),0.0,0.0};

            intervals[worstIdx] = left;
            intervals.insert(intervals.begin() + worstIdx + 1, right);

            compute_error(intervals[worstIdx]);
            compute_error(intervals[worstIdx+1]);

            // Recompute global metrics
            sum_base_T = std::complex<double>(0.,0.); sum_refined_T = std::complex<double>(0.,0.);
            sum_abs_diff_T = 0.0; worst_relerr_T = 0.0;
            for (auto& I : intervals) {
                sum_base_T += I.S_T; sum_refined_T += I.Sref_T;
                sum_abs_diff_T += std::abs(I.Sref_T - I.S_T);
                if (I.relerr > worst_relerr_T) worst_relerr_T = I.relerr;
            }
            denom_global_T = std::max(1.0, std::abs(sum_refined_T));
            global_abs_err_T = (sum_abs_diff_T / denom_global_T) / 15.0;
            splitsPerformedThetaB++;
        }

        // Accumulate with Richardson correction where refinement occurred
        std::complex<double> sumT(0.,0.), sumL(0.,0.); double sumT2=0.0, sumL2=0.0; int refinedAccepted=0;
        for (auto& I : intervals) {
            if (I.refined) {
                sumT  += I.Sref_T + (I.Sref_T - I.S_T)/15.0;
                sumT2 += I.Sref_T2 + (I.Sref_T2 - I.S_T2)/15.0;
                if (Qsqr > 0) {
                    sumL  += I.Sref_L + (I.Sref_L - I.S_L)/15.0;
                    sumL2 += I.Sref_L2 + (I.Sref_L2 - I.S_L2)/15.0;
                }
                refinedAccepted++;
            } else {
                sumT  += I.S_T; sumT2 += I.S_T2;
                if (Qsqr > 0) { sumL += I.S_L; sumL2 += I.S_L2; }
            }
        }
        out.F_T[ib] = sumT; out.F_T_sqr[ib] = sumT2;
        if (Qsqr > 0) { out.F_L[ib] = sumL; out.F_L_sqr[ib] = sumL2; }

        /*#pragma omp critical
        {
            std::cerr << "# theta_b evals for b=" << std::setprecision(6) << std::fixed << bval
                      << ": " << evals << " base " << baseIntervalsThetaB
                      << " intervals " << intervals.size() << " refined " << refinedAccepted
                      << " splits " << splitsPerformedThetaB
                      << " worst_relerr " << worst_relerr_T
                      << " global_abs_err " << global_abs_err_T
                      << " (cap " << maxEvalsCap << ")" << (evals >= maxEvalsCap ? " [CAP REACHED]" : "") << std::endl;
        }*/
    }

    // Integrate over b for total cross sections
    double sigma_T = 0.0, sigma_L = 0.0;
    for (int ib=0; ib<nbperp; ++ib) {
        double bval = out.b[ib];
        sigma_T += out.F_T_sqr[ib] * bval * db;
        if (Qsqr > 0) sigma_L += out.F_L_sqr[ib] * bval * db;
    }
    const double HBARC = 0.197327053; // GeV*fm
    sigma_T *= 1e7 * HBARC * HBARC / (8. * M_PI);
    sigma_L *= 1e7 * HBARC * HBARC / (8. * M_PI);
    out.sigma_T = sigma_T;
    out.sigma_L = sigma_L;
    return out;
}
