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
#include "subnucleon_config.hpp"
#include "nrqcd_wf.hpp"
#include <cmath>

using namespace std;

#include <complex>


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    num_of_averages = 1;
    zlimit=0.00000001;
    MAXR=10.*5.068;
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
    
    double Rg = 1.0; //std::pow(2.0, 2.0*lambda+3)/std::sqrt(M_PI) * gsl_sf_gamma(lambda+5.0/2.0)/gsl_sf_gamma(lambda+4.0);
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
    double mv; // mass of vector meson
    double daughter_mass;
    double root_snn;
    double mA; // target mass
    double BigP; //BigP = 0.5 * (p1-p2)
    double theta_BigP;
    double B;  // impact parameters between two nuclei
    double theta_B;
    double RA;
    int A;
    int Z;
    bool DacayToScalar;
    Polarization polarization;
    int M12reim;
};

double Integra_P_and_B( double *vec, size_t dim, void* par);
double Inthelperf_amplitude_and_F_mc( double *vec, size_t dim, void* par);
double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);
double Inthelperf_amplitude_mc_fixed_b( double *vec, size_t dim, void* par);
double Integra_P_and_B_mc( double *vec, size_t dim, void* par);

double Inthelperf_amplitude_soft_photon_mc( double *vec, size_t dim, void* par);

double Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, Polarization pol)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization=pol;

    
    
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
        lower[4]=zlimit; // Min z
        upper[4]=1.0 - lower[4];    // Max z
    }
    
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 10.*5.068 ; // Max b
    upper[1] = 2.0*M_PI;
    upper[2] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[3] = 2.0*M_PI;
    
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 4;
    if (!FACTORIZE_ZINT)
        F.dim = 5;
    F.params = &helper;
    
    double result,error;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        //cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        gsl_monte_vegas_free(s);
    }
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    delete lower;
    delete upper;
    
    return result;
    
}
double Diffraction::ScatteringAmplitude_at_fixed_b(double xpom, double Qsqr, double b, double theta_b, Polarization pol)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.polarization=pol;
    helper.b = b;
    helper.theta_b = theta_b;
    // Do MC integral over impact parameters and dipole sizes
    // Currently hardcoded parameters for jpsi and gold:
    // Impact parameter up to 100 GeV^-1
    // Dipole size up to 10 GeV^-1
    // MC integration parameters: b, theta_b, r, theta_r, Z
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
        lower[2]=zlimit; // Min z
        upper[2]=1.0 - lower[2];    // Max z
    }
    
    lower[0] = lower[1] = 0;
    upper[0] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[1] = 2.0*M_PI; // max theta_r
    
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc_fixed_b;
    F.dim = 2;
    if (!FACTORIZE_ZINT)
        F.dim = 3;
    F.params = &helper;
    
    double result,error;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        //cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        gsl_monte_vegas_free(s);
    }
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    delete lower;
    delete upper;
    
    return result;
    
}

double Diffraction::Soft_photon_ScatteringAmplitude(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, 
           double theta_BigP, int Z, double t, double RA, double Low, double High, double daughter_mass, bool DacayToScalar) {
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization=pol;
    
    helper.mv = mv; // the mass of vector meson
    helper.t = t;
    helper.root_snn = root_snn;
    helper.theta_BigP = theta_BigP;
    helper.Z = Z;
    helper.RA = RA;
    helper.daughter_mass = daughter_mass;
    helper.DacayToScalar = DacayToScalar;

    // Do MC integral over impact parameters and dipole sizes
    // Currently hardcoded parameters for jpsi and gold:
    // Impact parameter up to 100 GeV^-1
    // Dipole size up to 10 GeV^-1
    // MC integration parameters: b, theta_b( = theta_p), r, theta_r, Z, rT, theta_rT,  theta_B, B, Big_P
    double *lower, *upper;
    lower = new double[10];
    upper = new double[10];
    
    lower[0]=lower[1]=lower[2]=lower[3]=0.;
    lower[4] = zlimit;
    lower[5] = lower[6] = lower[7] = 0.0; 
    lower[8] = 2.* RA; // B

    upper[0] = 10.*5.068 ; // Max b
    upper[1] = 2.0*M_PI; //Max theta_b
    upper[2] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[3] = 2.0*M_PI;
    upper[4] = 1.0 - lower[4]; // MAX Z
    upper[5] = 10.*5.068; // max rT
    upper[6] = 2.0*M_PI;
    upper[7] = 2.0*M_PI;
    upper[8] = 10.*RA; // Max B

    double Q2Low = Low * mv * mv;
    double lowbigp = 0.5*std::sqrt((t*Q2Low+Q2Low*Q2Low)/(t+Q2Low - t*cos(theta_BigP)*cos(theta_BigP)));
    double Q2High = High * mv * mv;
    double highbigp = 0.5*std::sqrt((t*Q2High+Q2High*Q2High)/(t+Q2High - t*cos(theta_BigP)*cos(theta_BigP)));
    lower[9] = lowbigp;   //  0.7 M^2 < Q^2 < 1.3 M^2, Big_p
    upper[9] = highbigp;  //  0.7/4 M^2 < P^2 < 1.3/4 M^2

    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_soft_photon_mc;
    F.dim = 10;
    F.params = &helper;
    
    double result,error;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        gsl_monte_vegas_free(s);
    }
    
    delete lower;
    delete upper;
    
    return result;
    
}

double Diffraction::b_Inte(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, double B, double theta_B, 
                           int Z, double t, double BigP, double theta_BigP, int M12reim)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.mv = mv; // the mass of vector meson
    helper.t = t;
    helper.root_snn = root_snn;
    helper.B = B;       // impact parameters between two nuclei
    helper.theta_B = theta_B;
    helper.Z = Z;
    helper.M12reim = M12reim;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.polarization=pol;
    helper.BigP = BigP;
    helper.theta_BigP = theta_BigP;
    // Do MC integral over b, theta_b, r, theta_r, z
    // MC integration parameters: B
    double *lower, *upper;
    lower = new double[5];
    upper = new double[5];
    
    lower[0]=lower[1]=lower[2]=lower[3]=0.;
    lower[4] = zlimit;

    upper[0] = 10.*5.068 ; // Max b
    upper[1] = 2.0*M_PI; //Max theta_b
    upper[2] = MAXR; //MAXR;//20; //0.5*5.068;  // Max r
    upper[3] = 2.0*M_PI;
    upper[4] = 1.0 - lower[4]; // MAX Z

    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_and_F_mc;
    F.dim = 5;
    F.params = &helper;
    
    double result, error;
    const double VEGAS_RESULT_ACCURACY_TARGET=0.01;
    const int MAXITER_VEGAS=7;
    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# Vegas interation for B and P " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.3 or std::abs(error/result) > VEGAS_RESULT_ACCURACY_TARGET) and iter < MAXITER_VEGAS));
        gsl_monte_vegas_free(s);
    }
    
    delete lower;
    delete upper;
    
    return result;
}


double Diffraction::Relative_P_abd_B_Inte(double mv, double root_snn, double theta_BigP, int Z, double t, double RA, double Low, double High, bool DacayToScalar)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.mv = mv; // the mass of vector meson
    helper.t = t;
    helper.root_snn = root_snn;
    helper.theta_BigP = theta_BigP;
    helper.Z = Z;
    helper.RA = RA;
    helper.DacayToScalar = DacayToScalar;

    double BigP; //BigP = 0.5 * (p1-p2)
    double B;  // impact parameters between two nuclei
    // Do MC integral over impact parameters B and BigP
    // MC integration parameters: B
    double *lower, *upper;
    lower = new double[2];
    upper = new double[2];
    
    lower[0] = 2.0 * RA;
    upper[0] = 10.* RA; // Max B
    //lower[1] = 0.3535533905932738 * mv;  //  0.5 M^2 < Q^2 < 1.5 M^2
    //upper[1] = 0.6123724356957945 * mv;  //  0.5/4 M^2 < P^2 < 1.5/4 M^2
    //lower[1] = std::sqrt(Low)/2. * mv;  //  0.7 M^2 < Q^2 < 1.3 M^2
    //upper[1] = std::sqrt(High)/2. * mv;  //  0.7/4 M^2 < P^2 < 1.3/4 M^2
    double Q2Low = Low * mv * mv;
    double lowbigp = 0.5*std::sqrt((t*Q2Low+Q2Low*Q2Low)/(t+Q2Low - t*cos(theta_BigP)*cos(theta_BigP)));
    double Q2High = High * mv * mv;
    double highbigp = 0.5*std::sqrt((t*Q2High+Q2High*Q2High)/(t+Q2High - t*cos(theta_BigP)*cos(theta_BigP)));
    lower[1] = lowbigp;   //  0.7 M^2 < Q^2 < 1.3 M^2, Big_p
    upper[1] = highbigp;  //  0.7/4 M^2 < P^2 < 1.3/4 M^2
    
    gsl_monte_function F;
    F.f = &Integra_P_and_B;
    F.dim = 2;
    F.params = &helper;
    
    double result, error;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# Vegas interation for B and P " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        gsl_monte_vegas_free(s);
    }
    
    delete lower;
    delete upper;
    
    return result;
    
}

double Diffraction::Relative_P_abd_B_Inte_mc(double xpom, double Qsqr, Polarization pol, double mv, double root_snn, double theta_BigP, 
                                             int Z, double t, double RA, double Low, double High, bool DacayToScalar)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.mv = mv; // the mass of vector meson
    helper.t = t;
    helper.root_snn = root_snn;
    helper.theta_BigP = theta_BigP;
    helper.Z = Z;
    helper.RA = RA;
    helper.DacayToScalar = DacayToScalar;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization=pol;

    double BigP; //BigP = 0.5 * (p1-p2)
    double B;  // impact parameters between two nuclei
    // Do MC integral over impact parameters B, theta_B and BigP
    // MC integration parameters: B, theta_B, BigP
    double *lower, *upper;
    lower = new double[3];
    upper = new double[3];
    
    lower[0] = 2.0 * RA;
    upper[0] = 10.* RA; // Max B
    lower[1] = 0.0;
    upper[1] = 2.*M_PI; // Max theta_B
    double Q2Low = Low * mv * mv;
    double lowbigp = 0.5*std::sqrt((t*Q2Low+Q2Low*Q2Low)/(t+Q2Low - t*cos(theta_BigP)*cos(theta_BigP)));
    double Q2High = High * mv * mv;
    double highbigp = 0.5*std::sqrt((t*Q2High+Q2High*Q2High)/(t+Q2High - t*cos(theta_BigP)*cos(theta_BigP)));
    lower[2] = lowbigp;   //  0.7 M^2 < Q^2 < 1.3 M^2, Bigp
    upper[2] = highbigp;  //  0.7/4 M^2 < P^2 < 1.3/4 M^2
    
    gsl_monte_function F;
    F.f = &Integra_P_and_B_mc;
    F.dim = 3;
    F.params = &helper;
    
    double result, error;

       const double VEGAS_RESULT_ACCURACY_TARGET=0.01;
    const int MAXITER_VEGAS=7;
 
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        //cout << "# vegas warmup " << result << " +/- " << error << endl;
        int iter=0;
        do
        {
            iter++;
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# Vegas interation for B and P " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (iter < 2 or ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.3 or std::abs(error/result) > VEGAS_RESULT_ACCURACY_TARGET) and iter < MAXITER_VEGAS));
        gsl_monte_vegas_free(s);
    }
    
    delete lower;
    delete upper;
    
    return result;
    
}

double Inthelperf_amplitude_z(double z, void* p);

double Integra_P_and_B( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->B = vec[0];
    helper->BigP = vec[1];
    double alpha_em = 0.0073;
    double e_charge = 0.303;// GeV
    double BW_Gamma; // GeV, The Breit-Wigner width of the J/Psi From PDG. 
    double mproton = 0.938;// GeV
    double HbarC = 1.0; // GeV . fm
    double gamma = helper->root_snn/2./mproton;
    double Eq2phi02;
    double BW_prefactor;
    double bracket;
    double daughter_mass = 0.139;// GeV

    double nwB = helper->Z*helper->Z * alpha_em * helper->mv*helper->mv/4./M_PI/M_PI/gamma/gamma *  // w = mv/2*exp(-y), y=  0.0;
                 gsl_sf_bessel_Knu(1.0, helper->mv / 2. * helper->B  / gamma) *
                 gsl_sf_bessel_Knu(1.0, helper->mv / 2. * helper->B  / gamma);

    
    double N0tilde = helper->B * (1. - gsl_sf_bessel_J0(helper->B*std::sqrt(helper->t))) * nwB;
    double N2tilde = helper->B * gsl_sf_bessel_Jn(2.0, helper->B*std::sqrt(helper->t)) * nwB;
    if (helper->DacayToScalar) { // rho -> pi+ + pi-
        bracket = helper->BigP * helper->BigP / 2. * N0tilde + helper->BigP * helper->BigP/2. * N2tilde * cos( 2. * helper->theta_BigP);
        BW_prefactor = 12.24 * 12.24; // frhopipi = 12.24
        BW_Gamma = 0.156;//GeV, rho -> pipi GeV
        daughter_mass = 0.139;// GeV pion
    } else { // J/Psi -> mu+ + mu-
        BW_Gamma = 9.3e-05; // GeV, The Breit-Wigner width of the J/Psi From PDG. 
        Eq2phi02 = BW_Gamma * helper->mv * helper->mv /16./M_PI/alpha_em/alpha_em;
        BW_prefactor = 24. * pow(e_charge, 4) * Eq2phi02 / helper->mv;
        daughter_mass = 0.1056583745;// GeV muon 
        bracket = (1.-2.*helper->BigP * helper->BigP / helper->mv / helper->mv) * N0tilde - 
                   2.*helper->BigP * helper->BigP / helper->mv / helper->mv * N2tilde * cos( 2. * helper->theta_BigP); //theta_ p = 0.0
    }
    /*
    double Qsquare = 2. * (-0.25*helper->t + helper->BigP*helper->BigP + 
                     std::sqrt(pow(-0.5*sqrt(helper->t) + helper->BigP*cos(helper->theta_BigP), 2) +
                     pow(helper->BigP*sin(helper->theta_BigP), 2)) *
                     std::sqrt(pow(0.5*sqrt(helper->t) + helper->BigP*cos(helper->theta_BigP), 2) +
                     pow(helper->BigP*sin(helper->theta_BigP), 2)) 
                     );
    */
    double Qsquare = -0.5*helper->t + 2.*helper->BigP*helper->BigP + 2. *
                     std::sqrt( daughter_mass*daughter_mass + 0.25*helper->t + helper->BigP*helper->BigP -
                                std::sqrt(helper->t)*helper->BigP*cos(helper->theta_BigP) ) *
                     std::sqrt( daughter_mass*daughter_mass + 0.25*helper->t + helper->BigP*helper->BigP +
                                std::sqrt(helper->t)*helper->BigP*cos(helper->theta_BigP) );

    //double Qsquare = 4.* helper->BigP * helper->BigP;
    double Big_int = BW_prefactor/ pow(M_PI, 4)/ 16. * bracket / (helper->mv * helper->mv * BW_Gamma * BW_Gamma + 
                     (Qsquare - helper->mv * helper->mv) * (Qsquare - helper->mv * helper->mv));
    Big_int = std::max(Big_int, 0.0);
    return Big_int * helper->BigP / 2.; //helper->BigP / 2. is from Jacobians
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
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z, helper->polarization);
    
}

double Inthelperf_amplitude_mc_fixed_b( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->r = vec[0];
    helper->theta_r = vec[1];
    
    double z = 0.5;// Put z=0.5 as it sets b to the geometric average of quarks
    
    if (!FACTORIZE_ZINT)
        z = vec[2];
        
    return helper->diffraction->ScatteringAmplitudeIntegrand_fixed_b(helper->xpom, helper->Qsqr, helper->r, helper->theta_r, helper->b, helper->theta_b, z, helper->polarization);
}

double Inthelperf_amplitude_and_F_mc( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    double HbarC = 1.0; // GeV . fm
    helper->b = vec[0];
    helper->theta_b=vec[1];
    helper->r = vec[2];
    helper->theta_r = vec[3];
    helper->z = vec[4];
    //helper->z = 0.5;
    std::complex<double> imag(0.,1.);
    std::complex<double> Mx(1.,0.);
    std::complex<double> My(1.,0.);
    double alpha_em = 0.0073;
    double delta = std::sqrt(helper->t);
    double Bx = helper->B * cos(helper->theta_B);
    double By = helper->B * sin(helper->theta_B);
    double bx = helper->b * cos(helper->theta_b);
    double by = helper->b * sin(helper->theta_b);
    double Px = helper->BigP * cos(helper->theta_BigP);
    double Py = helper->BigP * sin(helper->theta_BigP);
    double mproton = 0.938;// GeV
    double gamma = helper->root_snn/2./mproton;
    double K_prefactor = helper->mv / 2.  / gamma;
    double Ftilde_part = helper->Z * std::sqrt(alpha_em) * helper->mv/2./M_PI /gamma; // w = mv/2*exp(-y), y=  0.0;
                    
    std::complex<double> niAtilde = helper->diffraction->ScatteringAmplitudeIntegrand_reim(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, helper->z, helper->polarization);
    std::complex<double> M0_part = niAtilde * Ftilde_part;
    double b_m_B = std::sqrt((Bx-bx)*(Bx-bx) + (By-by)*(By-by)) + 1.e-20;
    double b_p_B = std::sqrt((Bx+bx)*(Bx+bx) + (By+by)*(By+by)) + 1.e-20;
        Mx =  Bx * M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_m_B) / b_m_B - 
                                  M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_m_B) / b_m_B * bx - 
                                  (Bx * M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_p_B) / b_p_B + 
                                   M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_p_B) / b_p_B * bx) * 
                                  std::exp(-imag*(helper->B * delta * cos(helper->theta_B)));
        My =  By * M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_m_B) / b_m_B - 
                                  M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_m_B) / b_m_B * by - 
                                  (By * M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_p_B) / b_p_B + 
                                   M0_part * gsl_sf_bessel_Knu(1.0, K_prefactor * b_p_B) / b_p_B * by) * 
                                  std::exp(-imag*(helper->B * delta * cos(helper->theta_B)));
    //std::complex<double> MdotP = Mx * Px + My * Py;
    
    //return std::abs(MdotP) * std::abs(MdotP);
    double res = 0.0;
    if (helper->M12reim == 1) res = Mx.real() * Px + My.real() * Py;
    if (helper->M12reim == 2) res = Mx.imag() * Px + My.imag() * Py;
    return res;
}

double Integra_P_and_B_mc( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    helper->B = vec[0];
    helper->theta_B = vec[1];
    helper->BigP = vec[2];
    double alpha_em = 0.0073;
    double e_charge = 0.303;// GeV
    double BW_Gamma; // GeV, The Breit-Wigner width of the J/Psi From PDG. 
    double mproton = 0.938;// GeV
    double HbarC = 1.0; // GeV . fm
    double gamma = helper->root_snn/2./mproton;
    double Eq2phi02;
    double BW_prefactor;
    double bracket;
    double daughter_mass = 0.139;// GeV

    double Mreal = helper->diffraction->b_Inte(helper->xpom, helper->Qsqr, helper->polarization, helper->mv, helper->root_snn, 
                                             helper->B, helper->theta_B, helper->Z, helper->t, helper->BigP, helper->theta_BigP, 1);
    double Mimag = helper->diffraction->b_Inte(helper->xpom, helper->Qsqr, helper->polarization, helper->mv, helper->root_snn, 
                                             helper->B, helper->theta_B, helper->Z, helper->t, helper->BigP, helper->theta_BigP, 2);
    if (helper->DacayToScalar) { // rho -> pi+ + pi-
        bracket = Mreal * Mreal + Mimag * Mimag;
        BW_prefactor = 12.24 * 12.24; // frhopipi = 12.24
        BW_Gamma = 0.156;//GeV, rho -> pipi GeV
        daughter_mass = 0.139;// GeV pion
    } else { // J/Psi -> mu+ + mu-
        daughter_mass = 0.1056583745;// GeV muon 
        bracket = 0.0; //theta_ p = 0.0
    }
    double Qsquare = -0.5*helper->t + 2.*helper->BigP*helper->BigP + 2. *
                     std::sqrt( daughter_mass*daughter_mass + 0.25*helper->t + helper->BigP*helper->BigP -
                                std::sqrt(helper->t)*helper->BigP*cos(helper->theta_BigP) ) *
                     std::sqrt( daughter_mass*daughter_mass + 0.25*helper->t + helper->BigP*helper->BigP +
                                std::sqrt(helper->t)*helper->BigP*cos(helper->theta_BigP) );

    double Big_int = BW_prefactor/ pow(M_PI, 4)/ 16. * bracket / (helper->mv * helper->mv * BW_Gamma * BW_Gamma + 
                     (Qsquare - helper->mv * helper->mv) * (Qsquare - helper->mv * helper->mv));
    Big_int = std::max(Big_int, 0.0);
    //std::cout << " Big_int * helper->BigP / 2. = " << Big_int * helper->BigP / 2. <<std::endl;
    return Big_int * helper->BigP / 2.; //helper->BigP / 2. is from Jacobians
}

double Inthelperf_amplitude_soft_photon_mc( double *vec, size_t dim, void* par)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)par;
    double b = vec[0];
    double theta_b=vec[1];
    double r = vec[2];
    double theta_r = vec[3];
    double frac_z  = vec[4];
    double rT  = vec[5];
    double theta_rT  = vec[6];
    double theta_B  = vec[7];
    double B  = vec[8];
    double Big_P  = vec[9];
    
    // MC integration parameters: b, theta_b( = theta_p), r, theta_r, frac_z, rT, theta_rT,  theta_B, B, Big_P

    double alpha_em = 0.0073;
    double e_charge = 0.303;// GeV
    double BW_Gamma; // GeV, The Breit-Wigner width of the J/Psi From PDG. 
    double mproton = 0.938;// GeV
    double HbarC = 1.0; // GeV . fm
    double gamma = helper->root_snn/2./mproton;
    double Eq2phi02;
    double BW_prefactor;
    double bracket;

    auto niA_b      = helper->diffraction->ScatteringAmplitude_noexp_Integrand(helper->xpom, helper->Qsqr, helper->t, r, theta_r, b, theta_b, 0., 0., 0., 0., frac_z, helper->polarization);
    auto niA_bpB    = helper->diffraction->ScatteringAmplitude_noexp_Integrand(helper->xpom, helper->Qsqr, helper->t, r, theta_r, b, theta_b, B, theta_B, 0., 0., frac_z, helper->polarization);
    auto niA_bprT   = helper->diffraction->ScatteringAmplitude_noexp_Integrand(helper->xpom, helper->Qsqr, helper->t, r, theta_r, b, theta_b, rT, theta_rT, 0., 0., frac_z, helper->polarization);
    auto niA_bpBprT = helper->diffraction->ScatteringAmplitude_noexp_Integrand(helper->xpom, helper->Qsqr, helper->t, r, theta_r, b, theta_b, B, theta_B, rT, theta_rT, frac_z, helper->polarization);

    // Sudakov Factor
    double c2 = 2.*std::log(helper->mv/helper->daughter_mass) - 2.772588722239781;// -4log(2);
    double c4 = 2.*std::log(helper->mv/helper->daughter_mass) - 4.;
    double mur = 1.1229189671354916 / rT; //2.*std::exp(-0.5772156649)
    double Sudakov = (1.-2.*alpha_em*c2/M_PI*cos(2.*theta_rT - 2.*helper->theta_BigP) + 
                      alpha_em*c4/M_PI * cos(4.*theta_rT - 4.*helper->theta_BigP)) *
                      exp(-1.*alpha_em/M_PI * 2.*std::log(helper->mv/helper->daughter_mass) * 2.*std::log(Big_P/mur));

    double bracket2 = (niA_b.first - niA_bpB.first) * (niA_bprT.first - niA_bpBprT.first) 
                    + (niA_b.second - niA_bpB.second) * (niA_bprT.second - niA_bpBprT.second); // 0.006332573977646111 = 1/4/M_PI/M_PI

    double nwB = helper->Z*helper->Z * alpha_em * helper->mv*helper->mv/4./M_PI/M_PI/gamma/gamma *  // w = mv/2*exp(-y), y=  0.0;
                 gsl_sf_bessel_Knu(1.0, helper->mv / 2. * B  / gamma) *
                 gsl_sf_bessel_Knu(1.0, helper->mv / 2. * B  / gamma);

    if (helper->DacayToScalar) { // rho -> pi+ + pi-
        bracket = Big_P * Big_P / 2. * (cos(2.*theta_B - 2.*helper->theta_BigP)+1.);
        BW_prefactor = 12.24 * 12.24; // frhopipi = 12.24
        BW_Gamma = 0.156;//GeV, rho -> pipi GeV
    } else { // J/Psi -> mu+ + mu-
        BW_Gamma = 9.3e-05; // GeV, The Breit-Wigner width of the J/Psi From PDG. 
        Eq2phi02 = BW_Gamma * helper->mv * helper->mv /16./M_PI/alpha_em/alpha_em;
        BW_prefactor = 24. * pow(e_charge, 4) * Eq2phi02 / helper->mv;
        bracket = 1.-2.*Big_P*Big_P / helper->mv / helper->mv -2.*Big_P*Big_P / helper->mv / helper->mv *
                   cos(2.*theta_B - 2.*helper->theta_BigP); //theta_ p = 0.0
    }
    double Qsquare = 2. * (-0.25*helper->t + Big_P*Big_P + 
                     std::sqrt(pow(-0.5*sqrt(helper->t) + Big_P*cos(helper->theta_BigP), 2) +
                     pow(Big_P*sin(helper->theta_BigP), 2)) *
                     std::sqrt(pow(0.5*sqrt(helper->t) + helper->BigP*cos(helper->theta_BigP), 2) +
                     pow(Big_P*sin(helper->theta_BigP), 2)) 
                     );
    
    //double Qsquare = 4.* helper->BigP * helper->BigP;
    double Big_int = BW_prefactor/ pow(M_PI*2., 9)/ 2. / (helper->mv * helper->mv * BW_Gamma * BW_Gamma + 
                     (Qsquare - helper->mv * helper->mv) * (Qsquare - helper->mv * helper->mv)) *
                     nwB * bracket * Sudakov * cos(std::sqrt(helper->t) * rT  * cos(theta_rT)) * 
                     bracket2;
    Big_int = std::max(Big_int, 0.0);
    return Big_int * Big_P / 2.; //helper->BigP / 2. is from Jacobians
}


double Inthelperf_amplitude_z(double z, void* p)
{
    Inthelper_amplitude *helper = (Inthelper_amplitude*)p;
    
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z);
}

double Diffraction::ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol)
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
    z=0.5; // Do not use z when calcualting antiquark/quark positions, just b is geometric mean
    
    
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    
    z = tmpz;
    
    double delta = std::sqrt(t);
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    
    std::complex<double> result = 2.0*r*b; // r and b from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);
    
    if (FACTORIZE_ZINT)
    {
        if (wavef->WaveFunctionType() != "NRQCD")
        {
            //PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
            cerr << "FACTORIZE_ZINT currently only works with NRQCD wf" << endl;
            return 0;
        }
        // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
        if (pol == T)
            result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
        else
            result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
            
        result *= std::exp(-imag*(b*delta*std::cos(theta_b)))*amp;
        
    }
    else
    {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
        
        // This integrand is now not integrated over z
        std::complex<double> exponent = std::exp( -imag* ( b*delta*std::cos(theta_b) - (0.5 - z)*r*delta*std::cos(theta_r)  )  );
        
        result *= amp * exponent;
        
    }
    
	double res=0;
    if (REAL_PART)
        res = result.real();
    else
        res = result.imag();
    
    if (std::isnan(res) or std::isinf(res))
    {
        cerr << "Amplitude integral is " << res << " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
        exit(1);
    }
    
    return res;
}

double Diffraction::ScatteringAmplitudeIntegrand_fixed_b(double xpom, double Qsqr,  double r, double theta_r, double b, double theta_b, double z, Polarization pol)
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
    z=0.5; // Do not use z when calcualting antiquark/quark positions, just b is geometric mean
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    z = tmpz;
    double delta = 0.0;
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    std::complex<double> result = 2.0*r; // r from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);
    if (FACTORIZE_ZINT)
    {
        if (wavef->WaveFunctionType() != "NRQCD")
        {
            //PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
            cerr << "FACTORIZE_ZINT currently only works with NRQCD wf" << endl;
            return 0;
        }
        // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
        if (pol == T)
            result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
        else
            result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
        result *= std::exp(-imag*(b*delta*std::cos(theta_b)))*amp;
    }
    else
    {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
        // This integrand is now not integrated over z
        //std::complex<double> exponent = std::exp( -imag* ( b*delta*std::cos(theta_b) - (0.5 - z)*r*delta*std::cos(theta_r)  )  );
        result *= amp;
    }
	double res=0;
    if (REAL_PART)
        res = result.real();
    else
        res = result.imag();
    
    if (std::isnan(res) or std::isinf(res))
    {
        cerr << "Amplitude integral is " << res << " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr <<  " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
        exit(1);
    }
    return res;
}

std::complex<double> Diffraction::ScatteringAmplitudeIntegrand_reim(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol)
{ 
    
        
    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    
    // If do like Lappi, Mantysaari: set z=1/2 here
    
    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    double HbarC = 1.0; // GeV . fm
    // q and antiq positions
    double tmpz = z;
    z=0.5; // Do not use z when calcualting antiquark/quark positions, just b is geometric mean

    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    z = tmpz;
    double delta = std::sqrt(t);
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    
    std::complex<double> result = 2.0*r*b; // r and b from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);
    
    if (FACTORIZE_ZINT)
    {
        if (wavef->WaveFunctionType() != "NRQCD")
        {
            //PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
            cerr << "FACTORIZE_ZINT currently only works with NRQCD wf" << endl;
            return 0;
        }
        // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
        if (pol == T)
            result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
        else
            result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
            
        result *= std::exp(-imag*(b*delta*std::cos(theta_b)))*amp;
        
    }
    else
    {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
        
        // This integrand is now not integrated over z
        std::complex<double> exponent = std::exp( -imag* ( b*delta*std::cos(theta_b) - (0.5 - z)*r*delta*std::cos(theta_r) )  );
        
        result *= amp * exponent;
        
    }
    if (std::isnan(result.real()) or std::isinf(result.imag()))
    {
        cerr << "Amplitude integral is " << result.real() << "  " << result.imag() <<  " dipole " << amp << " xpom=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
        exit(1);
    }
    return result;
}

std::pair<double, double> Diffraction::ScatteringAmplitude_noexp_Integrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double bp, double theta_bp, double bp2, double theta_bp2, double z, Polarization pol)
{ 
    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    
    // If do like Lappi, Mantysaari: set z=1/2 here
    
    double bx = b*cos(theta_b) - bp*cos(theta_bp) - bp2*cos(theta_bp2);
    double by = b*sin(theta_b) - bp*sin(theta_bp) - bp2*sin(theta_bp2);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    
    // q and antiq positions
    double tmpz = z;
    z=0.5; // Do not use z when calcualting antiquark/quark positions, just b is geometric mean
    
    
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    
    z = tmpz;
    
    double delta = std::sqrt(t);
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    
    std::complex<double> result = 2.0*r*b; // r and b from Jacobians, 2 as we have written sigma_qq = 2 N
    std::complex<double> imag(0,1);
    
    if (FACTORIZE_ZINT)
    {
        if (wavef->WaveFunctionType() != "NRQCD")
        {
            //PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta)
            cerr << "FACTORIZE_ZINT currently only works with NRQCD wf" << endl;
            return std::make_pair(0.0, 0.0);;
        }
        // Note 1/(4pi) is included in the z integral measure in PsiSqr_T_intz
        if (pol == T)
            result *= ((NRQCD_WF*)wavef)->PsiSqr_T_intz(Qsqr, r, delta, theta_r);
        else
            result *= ((NRQCD_WF*)wavef)->PsiSqr_L_intz(Qsqr, r, delta,theta_r);
            
        //result *= std::exp(-imag*(b*delta*std::cos(theta_b)))*amp;
        result *= amp;
    }
    else
    {
        if (pol == T)
            result *= wavef->PsiSqr_T(Qsqr, r, z)/(4.0*M_PI); // Wavef
        else
            result *= wavef->PsiSqr_L(Qsqr, r, z)/(4.0*M_PI);
        
        // This integrand is now not integrated over z
        //std::complex<double> exponent = std::exp( -imag* ( b*delta*std::cos(theta_b) - (0.5 - z)*r*delta*std::cos(theta_r)  )  );
        
        result *= amp;// * exponent;
        
    }
    
	double res=0;
    /*
    if (REAL_PART)
        res = result.real();
    else
        res = result.imag();
    */
    if (std::isnan(result.real()) or std::isinf(result.imag()))
    {
        cerr << "Amplitude integral is " << result.real() << "  " << result.imag() <<  " dipole " << amp << " xp=" << xpom << " Q^2=" << Qsqr << " t="<< t << " r=" << r << " theta_r="<<theta_r << " b="<< b << "theta_b="<< theta_b << " z=" << z << endl;
        exit(1);
    }
    
    return std::make_pair(result.real(), result.imag());

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
    int status = gsl_integration_qag(&f, 0, 100, 0, 0.0001, INTPOINTS_ROTSYM, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#bint failed, result " << result << " relerror " << error << " t " <<t << endl;
    
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
    int status = gsl_integration_qag(&f, std::log(1e-6), std::log(50), 0, 0.0001, INTPOINTS_ROTSYM, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#Rint failed, result " << result << " relerror " << error << " b " << b << " t " << par->t << endl;
    
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
            cerr << "#zint failed, result " << result << " relerror " << error << " b " << par->b << " t " <<par->t << endl;
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
