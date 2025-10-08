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

        if (ShowVegasIterations())
            cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5,
                                      global_rng, s, &result, &error);
            if (ShowVegasIterations())
                cout << "# Vegas interation " << result << " +/- " << error
                     << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            if (iter>10)
                break;
        } while (iter < 2 or ( std::abs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 or std::abs(error/result) > MCINTACCURACY));
        gsl_monte_vegas_free(s);
    }
    delete lower;
    delete upper;
    return result;
}


double Inthelperf_amplitudeF_mc(double *vec, size_t dim, void* par);
double Inthelperf_amplitudeF_integrated_mc(double *vec, size_t dim, void* par);

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
    // MC integration parameters: r, theta_r, Z
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
                cout << "# Vegas interation " << result << " +/- " << error
                     << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
            if (iter>10)
                break;
        } while (iter < 2 or ( std::abs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 or std::abs(error/result) > MCINTACCURACY));
        gsl_monte_vegas_free(s);
    } else if (MCINT == ADAPTIVE_QUAD) {
        // Use trapezoidal rule with smart grid spacing
        int nr_points = 64;   // Much fewer points than MC for speed
        int ntheta_points = 32;  // Sufficient for angular integration
        int nz_points = (!FACTORIZE_ZINT) ? 32 : 1;

        result = 0.0;
        
        // r integration with logarithmic spacing
        double r_min = 1e-6;  // Avoid r=0 exactly  
        double r_max = upper[0];
        double log_r_min = log(r_min);
        double log_r_max = log(r_max);
        double dlog_r = (log_r_max - log_r_min) / nr_points;

        for (int i = 0; i < nr_points; i++) {
            // Logarithmic spacing with midpoint rule
            double log_r = log_r_min + (i + 0.5) * dlog_r;
            double r = exp(log_r);
            helper.r = r;
            
            // theta_r integration
            double dtheta = 2.0 * M_PI / ntheta_points;
            double theta_sum = 0.0;
            
            for (int j = 0; j < ntheta_points; j++) {
                double theta_r = (j + 0.5) * dtheta;
                helper.theta_r = theta_r;
                
                if (FACTORIZE_ZINT) {
                    double z = 0.5;
                    double integrand = helper.diffraction->ScatteringAmplitudeIntegrand2(
                        helper.xpom, helper.Qsqr, helper.r, helper.theta_r,
                        helper.b, helper.theta_b, z, helper.polarization, helper.real_part);
                    theta_sum += integrand;
                } else {
                    double zlimit = lower[2];
                    double z_min = zlimit;
                    double z_max = upper[2];
                    double dz = (z_max - z_min) / nz_points;
                    double z_sum = 0.0;
                    
                    for (int k = 0; k < nz_points; k++) {
                        double z = z_min + (k + 0.5) * dz;
                        double integrand = helper.diffraction->ScatteringAmplitudeIntegrand2(
                            helper.xpom, helper.Qsqr, helper.r, helper.theta_r,
                            helper.b, helper.theta_b, z, helper.polarization, helper.real_part);
                        z_sum += integrand;
                    }
                    theta_sum += z_sum * dz;
                }
            }
            // Include correct Jacobians: only r from log integration (cylindrical r already in integrand)
            result += theta_sum * r * dlog_r * dtheta;
        }
    }
    delete lower;
    delete upper;
    return result;
}

// Compute integral over theta_b of F(b,theta_b) 
double Diffraction::ScatteringAmplitudeFIntegrated(
    double xpom, double Qsqr, double b, Polarization pol, bool real_part) {

    if (MCINT == ADAPTIVE_QUAD) {
        int ntheta_b_points = 32;  // theta_b integration points
        double dtheta_b = 2.0 * M_PI / ntheta_b_points;
        double result = 0.0;

        for (int ib = 0; ib < ntheta_b_points; ib++) {
            double theta_b = (ib + 0.5) * dtheta_b;
            
            int nr_points = 64;
            int ntheta_points = 32;
            int nz_points = (!FACTORIZE_ZINT) ? 32 : 1;

            double f_contribution = 0.0;

            // r integration with logarithmic spacing
            double r_min = 1e-10;  // Avoid r=0 exactly
            double r_max = MAXR;
            double log_r_min = log(r_min);
            double log_r_max = log(r_max);
            double dlog_r = (log_r_max - log_r_min) / nr_points;

            for (int i = 0; i < nr_points; i++) {
                // Logarithmic spacing with midpoint rule
                double log_r = log_r_min + (i + 0.5) * dlog_r;
                double r = exp(log_r);

                // theta_r integration
                double dtheta = 2.0 * M_PI / ntheta_points;
                double theta_sum = 0.0;

                for (int j = 0; j < ntheta_points; j++) {
                    double theta_r = (j + 0.5) * dtheta;

                    if (FACTORIZE_ZINT) {
                        double z = 0.5;
                        double integrand = ScatteringAmplitudeIntegrand2(
                            xpom, Qsqr, r, theta_r, b, theta_b, z, pol, real_part);
                        theta_sum += integrand;
                    } else {
                        double zlimit = 0.00000001;
                        double z_min = zlimit;
                        double z_max = 1.0 - zlimit;
                        double dz = (z_max - z_min) / nz_points;
                        double z_sum = 0.0;
                        
                        for (int k = 0; k < nz_points; k++) {
                            double z = z_min + (k + 0.5) * dz;
                            double integrand = ScatteringAmplitudeIntegrand2(
                                xpom, Qsqr, r, theta_r, b, theta_b, z, pol, real_part);
                            z_sum += integrand;
                        }
                        theta_sum += z_sum * dz;
                    }
                }
                // Include correct Jacobians: only r from log integration (cylindrical r already in integrand)
                f_contribution += theta_sum * r * dlog_r * dtheta;
            }
            result += f_contribution * dtheta_b;
        }
        return result;
    } else {
        // For Monte Carlo methods, integrate over theta_b, r, theta_r, (z) directly
        Inthelper_amplitude helper;
        helper.diffraction = this;
        helper.xpom = xpom;
        helper.Qsqr = Qsqr;
        helper.t = 0.0;
        helper.b = b;
        helper.polarization = pol;
        helper.real_part = real_part;

        double *lower, *upper;
        int dim;

        if (FACTORIZE_ZINT) {
            dim = 3;  // theta_b, r, theta_r
            lower = new double[3];
            upper = new double[3];
        } else {
            dim = 4;  // theta_b, r, theta_r, z
            lower = new double[4];
            upper = new double[4];
            lower[3] = zlimit;
            upper[3] = 1.0 - lower[3];
        }

        lower[0] = 0; upper[0] = 2.0*M_PI;  // theta_b
        lower[1] = 0; upper[1] = MAXR;      // r  
        lower[2] = 0; upper[2] = 2.0*M_PI;  // theta_r

        gsl_monte_function F;
        F.f = &Inthelperf_amplitudeF_integrated_mc;
        F.dim = dim;
        F.params = &helper;

        double result, error;

        if (MCINT == VEGAS) {
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50,
                                      global_rng, s, &result, &error);

            int iter = 0;
            do {
                iter++;
                gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5,
                                          global_rng, s, &result, &error);
                if (iter > 10) break;
            } while (iter < 2 or (std::abs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 or std::abs(error/result) > MCINTACCURACY));

            gsl_monte_vegas_free(s);
        } else {
            gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
            gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS,
                                      global_rng, s, &result, &error);
            gsl_monte_miser_free(s);
        }
        delete[] lower;
        delete[] upper;
        return result;
    }
}

// Compute integral over theta_b of |F(b,theta_b)|^2
double Diffraction::ScatteringAmplitudeFSquaredIntegrated(
    double xpom, double Qsqr, double b, Polarization pol) {

    if (MCINT == ADAPTIVE_QUAD) {
        // Compute |F|^2 by integrating |F(theta_b)|^2 over theta_b
        // where F(theta_b) is computed using the same method as ScatteringAmplitudeF
        int ntheta_b_points = 32;
        double dtheta_b = 2.0 * M_PI / ntheta_b_points;
        double result = 0.0;

        for (int ib = 0; ib < ntheta_b_points; ib++) {
            double theta_b = (ib + 0.5) * dtheta_b;

            // Compute F_real and F_imag at this theta_b using the same grid as ScatteringAmplitudeF
            int nr_points = 64;
            int ntheta_points = 32;
            int nz_points = (!FACTORIZE_ZINT) ? 32 : 1;

            // Use logarithmic spacing for r
            double r_min = 1e-10;
            double r_max = MAXR;
            double log_r_min = log(r_min);
            double log_r_max = log(r_max);
            double dlog_r = (log_r_max - log_r_min) / nr_points;
            double dtheta = 2.0 * M_PI / ntheta_points;

            double f_real = 0.0, f_imag = 0.0;

            for (int i = 0; i < nr_points; i++) {
                double log_r = log_r_min + (i + 0.5) * dlog_r;
                double r = exp(log_r);

                double theta_sum_real = 0.0, theta_sum_imag = 0.0;
                for (int j = 0; j < ntheta_points; j++) {
                    double theta_r = (j + 0.5) * dtheta;

                    if (FACTORIZE_ZINT) {
                        double z = 0.5;
                        double integrand_real = ScatteringAmplitudeIntegrand2(
                            xpom, Qsqr, r, theta_r, b, theta_b, z, pol, true);
                        double integrand_imag = ScatteringAmplitudeIntegrand2(
                            xpom, Qsqr, r, theta_r, b, theta_b, z, pol, false);
                        theta_sum_real += integrand_real;
                        theta_sum_imag += integrand_imag;
                    } else {
                        double zlimit = 0.00000001;
                        double z_min = zlimit;
                        double z_max = 1.0 - zlimit;
                        double dz = (z_max - z_min) / nz_points;
                        double z_sum_real = 0.0, z_sum_imag = 0.0;

                        for (int k = 0; k < nz_points; k++) {
                            double z = z_min + (k + 0.5) * dz;
                            double integrand_real = ScatteringAmplitudeIntegrand2(
                                xpom, Qsqr, r, theta_r, b, theta_b, z, pol, true);
                            double integrand_imag = ScatteringAmplitudeIntegrand2(
                                xpom, Qsqr, r, theta_r, b, theta_b, z, pol, false);
                            z_sum_real += integrand_real;
                            z_sum_imag += integrand_imag;
                        }
                        theta_sum_real += z_sum_real * dz;
                        theta_sum_imag += z_sum_imag * dz;
                    }
                }
                f_real += theta_sum_real * r * dlog_r * dtheta;
                f_imag += theta_sum_imag * r * dlog_r * dtheta;
            }
            // Compute |F|^2 = F_real^2 + F_imag^2 at this theta_b
            double f_squared = f_real * f_real + f_imag * f_imag;
            result += f_squared * dtheta_b;
        }
        return result;
    } else {
        // For Monte Carlo, use the linear approach
        double real_part = ScatteringAmplitudeFIntegrated(xpom, Qsqr, b, pol, true);
        double imag_part = ScatteringAmplitudeFIntegrated(xpom, Qsqr, b, pol, false);
        return real_part * real_part + imag_part * imag_part;
    }
}

double Inthelperf_amplitude_z(double z, void* p);

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
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z, helper->polarization, helper->real_part);
    
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

double Inthelperf_amplitudeF_integrated_mc(double *vec, size_t dim, void* par) {
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->theta_b = vec[0];
    helper->r = vec[1];
    helper->theta_r = vec[2];

    double z = 0.5;// Put z=0.5 as it sets b to the geometric average of quarks

    if (!FACTORIZE_ZINT)
        z = vec[3];

    return helper->diffraction->ScatteringAmplitudeIntegrand2(
        helper->xpom, helper->Qsqr, helper->r, helper->theta_r,
        helper->b, helper->theta_b, z, helper->polarization, helper->real_part);
}

double Inthelperf_amplitude_z(double z, void* p)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)p;
    
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z);
}

double Diffraction::ScatteringAmplitudeIntegrand(
    double xpom, double Qsqr, double t, double r, double theta_r,
    double b, double theta_b, double z, Polarization pol, bool real_part) {

    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    // If do like Lappi, Mantysaari: set z=1/2 here

    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);

    // q and antiq positions
    // Note my convention is that b is the center of the dipole (geometric center), not center of mass (z weighted)
    // Consequently I get (0.5-z)r.Delta phase    
    double qx = bx + 0.5*rx; double qy = by + 0.5*ry;
    double qbarx = bx - 0.5*rx; double qbary = by - 0.5*ry;

    double delta = std::sqrt(t);

    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    std::complex<double> amp = dipole->ComplexAmplitude(xpom, x1, x2); //(amp_real, amp_imag);
    std::complex<double> result = 2.0*r*b; // r and b from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);

    if (FACTORIZE_ZINT) {
        if (wavef->WaveFunctionType() == "NRQCD") {
            // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
            if (pol == T)
                result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
            else
                result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
        } else {
            if (pol == T)
                result *= wavef->PsiSqr_T_intz(Qsqr, r);
            else
                result *= wavef->PsiSqr_L_intz(Qsqr, r);
        }
        if (delta > 0) {
            result *= std::exp(-imag*(b*delta*std::cos(theta_b)))*amp;
        } else {
            result *= amp;
        }
    } else {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);

        // This integrand is now not integrated over z
        std::complex<double> exponent(1, 0);
        if (delta > 0) {
            exponent = std::exp( -imag* ( b*delta*std::cos(theta_b) - (0.5 - z)*r*delta*std::cos(theta_r)  )  );
        }

        result *= amp * exponent;

    }

    double res = 0;
    // Note: As I'm using GSL to integrate and I use a routine that assumes a scalar function, this is quite inefficeint as 
    // everything above is computed separately for the real and imaginary parts. Performance could be optimized by using an 
    // integration routine supportin vector valued functions
    if (real_part)
        res = result.real();
    else
        res = result.imag();

    if (std::isnan(res) or std::isinf(res)) {
        cerr << "Amplitude integral is " << res << " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
        exit(1);
    }

    return res;
}


double Diffraction::ScatteringAmplitudeIntegrand2(
    double xpom, double Qsqr, double r, double theta_r,
    double b, double theta_b, double z, Polarization pol, bool real_part) { 

    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    // If do like Lappi, Mantysaari: set z=1/2 here

    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);

    // q and antiq positions
    double qx = bx + (1. - z) * rx;
    double qy = by + (1. - z) * ry;
    double qbarx = bx - z * rx;
    double qbary = by - z * ry;

    double x1[2] = {qx, qy};
    double x2[2] = {qbarx, qbary};
    std::complex<double> amp = dipole->ComplexAmplitude(xpom, x1, x2); //(amp_real, amp_imag);

    std::complex<double> result = 2.0*r; // r from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);

    if (FACTORIZE_ZINT) {
        if (wavef->WaveFunctionType() == "NRQCD") {
            double delta = 0.;
            // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
            if (pol == T)
                result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
            else
                result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
        } else {
            if (pol == T)
                result *= wavef->PsiSqr_T_intz(Qsqr, r);
            else
                result *= wavef->PsiSqr_L_intz(Qsqr, r);
        }
        result *= amp;
    } else {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
        result *= amp;
    }

    double res = 0;
    // Note: As I'm using GSL to integrate and I use a routine that assumes a scalar function, this is quite inefficeint as 
    // everything above is computed separately for the real and imaginary parts. Performance could be optimized by using an 
    // integration routine supportin vector valued functions
    if (real_part)
        res = result.real();
    else
        res = result.imag();
    return res;
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
