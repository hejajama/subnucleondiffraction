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
#include <fstream>
#include <sstream>
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
    
    InitializeIsInterpolator("./photon_kT_Isfun_LHC");
    

    
}


Diffraction::~Diffraction()
{
    delete Is_interpolator;
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
    bool xcomp; // If true, compute Mx ,otherwise My
    double B;
    double theta_B;
    Polarization polarization;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

double* Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, double B, double theta_B, Polarization pol)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization=pol;
    helper.xcomp = true;

    
    
    // b, theta_b, r, theta_r, z
    double *lower, *upper;
    lower = new double[5];
    upper = new double[5];
    lower[4]=zlimit; // Min z
    upper[4]=1.0 - lower[4];    // Max z
    
    
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 25*5.068; //MAXR;//20; //0.5*5.068;  // Max r
    upper[1] = 2.0*M_PI;
    upper[2]=MAXR;
    upper[3] = 2.0*M_PI;
    
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 5;

    F.params = &helper;
    
    helper.B=B;
    helper.theta_B=theta_B;



    // X comp
    helper.xcomp=true;
    
    double result_x,error_x;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result_x, &error_x);
        //cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result_x, &error_x);
        cout << "# vegas warmup " << result_x << " +/- " << error_x << endl;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result_x, &error_x);
            cout << "# Vegas interation " << result_x << " +/- " << error_x << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 and result_x != 0);
        gsl_monte_vegas_free(s);
    }
    
    // Y comp
    helper.xcomp=false;
    
    double result_y,error_y;

    
    if (MCINT == MISER)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(F.dim);
        gsl_monte_miser_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result_y, &error_y);
        //cout << "# Miser result " << result << " err " << error << " relerr " << std::abs(error/result) << endl;
        gsl_monte_miser_free(s);
    }
    else if (MCINT == VEGAS)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(F.dim);
        gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/50, global_rng, s, &result_y, &error_y);
        cout << "# vegas warmup " << result_y << " +/- " << error_y << endl;
        do
        {
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result_y, &error_y);
            cout << "# Vegas interation " << result_y << " +/- " << error_y << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
        } while (fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.5 and result_y != 0);
        gsl_monte_vegas_free(s);
    }
    
   
    
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    delete lower;
    delete upper;
    
    double *res = new double[2]; res[0]=result_x; res[1] = result_y;
    return res;
    
}



double Inthelperf_amplitude_z(double z, void* p);

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* p)
{
    // Integration variables: (choose momentum transfer || x axis)
    // B, theta_B, b, theta_b, r, theta_r, z
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    
    if (par->polarization != T)
    {
        cerr << "Only T supported atm..." << endl;
        exit(1);
    }
    
    Vec b (vec[0]*std::cos(vec[1]), vec[0]*std::sin(vec[1]));
    Vec r (vec[2]*std::cos(vec[3]), vec[2]*std::sin(vec[3]));
    double qt = std::sqrt(par->t);
    Vec q(qt, 0); // Choose along x axis
    double z = vec[4];
    
    
    
    
    double amp_real = par->diffraction->GetDipole()->Amplitude(par->xpom, b + r*0.5,  b - r*0.5);
    //double amp_imag = dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    double amp_imag=0; //todo...
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    
    // vec[2] = |r| jacobian
    complex<double> subamp = vec[2]*par->diffraction->GetWaveFunction()->PsiSqr_T(par->Qsqr, vec[2], z)/(4.0*M_PI) * amp;
    
    complex<double> result;
    complex<double> imag(0,1);
    
    // vec[0]=|b| is Jacobian
    double blen = vec[0]; double B = par->B; double th_B = par->theta_B;
    
    double cos_th_b = std::cos(vec[1]);
    double cos_th_B = std::cos(th_B);
    double sin_th_b = std::sin(vec[1]);
    double sin_th_B = std::sin(th_B);
    
    Vec Bv(B * cos_th_B, B*sin_th_B );
    
    result = vec[0] * std::exp(-imag*(q*b)) * subamp;
    

    if (par->xcomp == true)
    {
        result *= (blen*cos_th_b - B*cos_th_B)/( (b-Bv).Len() ) * par->diffraction->GetIsInterpolator()->Evaluate((b-Bv).Len())
        + (blen*cos_th_b + B*cos_th_B)/( (b+Bv).Len() ) * par->diffraction->GetIsInterpolator()->Evaluate((b+Bv).Len()) * std::exp(-imag*(q*Bv));
    }
    else
    {
        result *= (blen*sin_th_b - B*sin_th_B)/( (b-Bv).Len() ) * par->diffraction->GetIsInterpolator()->Evaluate((b-Bv).Len())
        + (blen*sin_th_b + B*sin_th_B)/( (b+Bv).Len() ) * par->diffraction->GetIsInterpolator()->Evaluate((b+Bv).Len()) * std::exp(-imag*(q*Bv));
    }
    
    double res;
    if (REAL_PART)
        res = result.real();
    else
        res = result.imag();
    
    return res;
}



/*
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


/// Interpolated is
void Diffraction::InitializeIsInterpolator(std::string datafile)
{
    vector<double> x;
    vector<double> y;
    
    cout << "# Reading file " << datafile << endl;
    
    std::ifstream infile(datafile);
    
    for (std::string line; getline(infile, line); )
    {
        std::stringstream ss(line);
        double a,b;
        ss >> a;
        ss >> b;
        x.push_back(a);
        y.push_back(b);
    }
    
    Is_interpolator = new Interpolator(x,y);
    Is_interpolator->SetFreeze(true);
    Is_interpolator->SetUnderflow(0);
    Is_interpolator->SetOverflow(0);
    
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
