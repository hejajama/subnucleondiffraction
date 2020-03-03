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

using namespace std;

#include <complex>


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    num_of_averages = 1;
    zlimit=0.00000001;
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

    
    Polarization polarization;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

bool DIJET = true;


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
    
    lower = new double[4];
    upper = new double[4];
   
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 5*5.068;  ; //1*5.068; //100; // Max b
    upper[1] = 2.0*M_PI;
    upper[2] = 5*5.068;//20; //0.5*5.068;  // Max r
    upper[3] = 2.0*M_PI;
    
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 4;
    F.params = &helper;
    
    double result,error;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
        cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result, &error);
        cout << "# vegas warmup " << result << " +/- " << error << endl;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            cout << "# Vegas interation " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 and result != 0);
        gsl_monte_vegas_free(s);
    }
    
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
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z, helper->polarization);
    
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


double Diffraction::ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double t, double r, double theta_r, double b, double theta_b, double z, Polarization pol)
{ 
    
        
    // Recall quark and gluon positions:
    // Quark: b + 1/2 r
    // Antiquark: b - 1/2 r

    
    double bx = b*cos(theta_b);
    double by = b*sin(theta_b);
    double rx = r*cos(theta_r);
    double ry = r*sin(theta_r);
    

    
    
    double qx = bx + rx/2.0; double qy = by + ry/2.0;
    double qbarx = bx - rx/2.0; double qbary = by - ry/2.0;

    
    double delta = std::sqrt(t);
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    
    
    if (DIJET)
    {
        double dphi = t;
        
        // Setup: select kinematics (jet rapidities = z1,z2)
        double z0=0.5; double z1=0.5;
        // Quark flavor, note that currently only one mass for quarks is supported
//        double mq=1.28;
        double mq=0.14;
        // sum_f e_f^2
        double efsum = std::pow(2.0/3.0,2.0); //charm
        double eps = std::sqrt(Qsqr*z0*z1+mq*mq);
        
        // Set jet pt's, vary angle
        /*
        double pt0  =0.5;
        double pt1 = 0.5;
        
        
        // Construct dot product
        // Note that now all angles are measured w.r.t. p0, which is set to point along the x axis
        // b . (p0 + p1)
        double b_dot_p0plusp1 = b*pt0*cos(theta_b) + b*pt1*cos(dphi - theta_b);
        // r . (p0 - p1)
        double r_dot_p0minusp1 =r*pt0*cos(theta_r) - r*pt1*cos(dphi - theta_r);
         */
        
        // Set Delta = \vec p1 + \vec p2, P = 1/2 (\vec p1 - \vec p2)
        //double Delta = 0.4; // Members of the class
        //double P = 2.2;     // k in 1511.07452
        
        // t is angle between P and delta
        // Now all angles are measured wrt delta, which we choose to point along the x axis
        double b_dot_p0plusp1 = b*Delta*std::cos(theta_b);
        double r_dot_p0minusp1 = r*P*std::cos(dphi - theta_r);
        
        
        std::complex<double> imag(0,1);
        std::complex<double> exponent = std::exp( -imag* ( b_dot_p0plusp1 + r_dot_p0minusp1/2.0 )  );
        
        complex<double> result;
        
        // NOTE
        // Here the full photon wave function coeffs are not impelmented.
        // For simple testing, we just take the part of the wave function that depends on
        // dipole size r
        
        // L
        if (pol == L)
        {
            result = r*b*exponent * amp * gsl_sf_bessel_K0(eps*r);
            // Longitudinal photon wave function = K_0(eps*r)
            
            // Constants, following (20) in 1511.07452
            // Not including (2pi)^2 delta(z_0+z_1 -1)
            const int Nc=3;
            const double alpha_em = 1.0/137.0;
            result *=Nc*alpha_em * efsum * z0*z1 * 8.0*z0*z1*eps*eps / std::pow(2.0, 8.0);
        }
        
        // T
        if (pol == T)
        {
            // Wave function is K_1(eps*r), or in mass term K_0(eps*r)
            // As transverse xs is \sim | polarizationvector*r |^2, so amplitude
            // is polarizatoinvector * r, and eventually we polvec.r polvec.r' = r.r'
            // So the amplitude is a vector, and we need to compute separately its x and y
            // components
            /// TODO: chec that factors of r are correct here!
            if (dijet_component == X)
                result = r*b*exponent * amp * eps*gsl_sf_bessel_K1(eps*r) * r*sin(theta_r) / r;
            else if (dijet_component == Y)
                result = r*b*exponent * amp * eps*r*gsl_sf_bessel_K1(eps*r) * r*cos(theta_r)/r;
            else // mass term
                result = mq*r*b*exponent * amp * gsl_sf_bessel_K0(eps*r);
            // Coeffs missing
            
        }
        
    // Note Jacobian r*b is included already
        
        if (REAL_PART)
            return result.real();
        else
            return result.imag();
        
        
    }
    
    cerr << "This code only calculates DIJET!" << endl;
    exit(1);

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
