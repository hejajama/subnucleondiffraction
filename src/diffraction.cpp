/*
 * Diffraction at sub-nucleon scale
 * Calculate diffractive cross sections
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#include "diffraction.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_gamma.h>
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
 *
 * This code only returns Im A, later used should average and take the square
 */
double Diffraction::CoherentCrossSection(double xpom, double Qsqr, double t)
{
    if (!CORRECTIONS)
        return ScatteringAmplitude(xpom, Qsqr, t);

    cerr << "Corrctions not implemented!" << endl;
    return 0;
    // Calculate also real part and skewedness
    // Very crude lowest order approximation for the derivative
    double amp = ScatteringAmplitude(xpom, Qsqr, t);
    double y = std::log(1.0/xpom);
    double yp = y + DELTA_Y;
    double amp2 = ScatteringAmplitude(std::exp(-yp), Qsqr, t);
    
    double lambda = std::log(amp2/amp)/DELTA_Y;
    // Skewedness, see e.g. 0712.2670
    double rg = std::pow( 2.0, 2.0*lambda+3.0)/std::sqrt(M_PI);
    rg *= gsl_sf_gamma(lambda + 5.0/2.0)/gsl_sf_gamma(lambda+4.0);
    
    // Real part
    beta = std::tan( M_PI * lambda / 2.0);
    
    return amp;
}

double Diffraction::GetBeta()
{
    return beta;
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
        lower[4]=0.000001; // Min z
        upper[4]=1.0 - lower[4];    // Max z
    }
    
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 100; // Max b
    upper[1] = 2.0*M_PI;
    upper[2] = 10;  // Max r
    upper[3] = 2.0*M_PI;
     
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 4;
    if (!FACTORIZE_ZINT)
        F.dim = 5;
    F.params = &helper;
    
    double result,error;
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, r, s, &result, &error);
        cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, r, s, &result, &error);
        cout << "# vegas warmup " << result << " +/- " << error << endl;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, r, s, &result, &error);
            cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 and result != 0);
        gsl_monte_vegas_free(s);
    }
    
    gsl_rng_free(r);
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    delete lower;
    delete upper;
    
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
    
    double z = 0.5;// Put z=0.5 as it sets b to the geometric average of quarks
    
    if (!FACTORIZE_ZINT)
        z = vec[4];
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z);
    
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
    
    if (Qsqr>0)
        cerr <<"Check that Q^2>0 works (different polarizations)s!" << endl;
    
    // Recall quark and gluon positions:
    // Quark: b + zr
    // Antiquark: b - (1-z) r
    
    // If do like Lappi, Mantysaari: set z=1/2 here
    
    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    
    // q and antiq positions
    //double tmpz = z;
    if (FACTORIZE_ZINT)
        z=0.5;      // Use b as geometric average, decouple zintegral
    //z=0.5;
    
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;

    //z = tmpz;
    double res = 0;
    
    res = 2.0 * r * b;
    
    double delta = std::sqrt(t);
    
    if (FACTORIZE_ZINT)
    {
        res *= wavef->PsiSqr_T_intz(Qsqr, r);   // Note: 4pi factor is in PsiSqr_T_intz function!
        
        if (REAL_PART)
            res *= std::cos( b*delta*std::cos(theta_b));    // Neglect z now
        else
            res *= -std::sin( b*delta*std::cos(theta_b));
    }
    else
    {
        res *= wavef->PsiSqr_tot(Qsqr, r, z)/(4.0*M_PI); // Wavef
        // As this integrand is now not integrated over z
        if (REAL_PART)
            res *=std::cos( b*delta*std::cos(theta_b) - (1.0 - z)*r*delta*std::cos(theta_r));
        else
            res *=-std::sin( b*delta*std::cos(theta_b) - (1.0 - z)*r*delta*std::cos(theta_r));
    }
    
    
    
    // Real part cos, imaginary part -sin
    //res *= std::cos( b*delta*std::cos(theta_b) - (1.0 - z)*r*delta*std::cos(theta_r));

    
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