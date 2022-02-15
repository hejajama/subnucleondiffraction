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
    
    if (KINEMATICS == LHC)
        InitializeIsInterpolator("./photon_kT_Isfun_LHC");
    else if (KINEMATICS == RHIC)
        InitializeIsInterpolator("./photon_kT_Isfun_RHIC");
    

    
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
    COMPONENT component;
    double B;
    double theta_B;
    bool real_part;
    Polarization polarization;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

std::vector<double> Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, double B, double theta_B, bool real_part, std::vector<COMPONENT> complist, Polarization pol)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.polarization=pol;
    helper.real_part = real_part;
    helper.B=B;
    helper.theta_B=theta_B;


    
    
    // b, theta_b, r, theta_r, z
    double *lower, *upper;
    lower = new double[5];
    upper = new double[5];
    lower[4]=zlimit; // Min z
    upper[4]=1.0 - lower[4];    // Max z
    
    
    lower[0]=lower[1]=lower[2]=lower[3]=0;
    upper[0] = 30*5.068; // max b
    upper[1] = 2.0*M_PI;
    upper[2]=MAXR;
    upper[3] = 2.0*M_PI;
    
    gsl_monte_function F;
    F.f = &Inthelperf_amplitude_mc;
    F.dim = 5;

    F.params = &helper;
    
    
    
    const double VEGAS_RESULT_ACCURACY_TARGET=0.01;
    const int MAXITER_VEGAS=7;
    
    std::vector<double> results;
    for (int i=0; i<complist.size(); i++)
    {
        helper.component = complist[i];
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
            gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS/5, global_rng, s, &result, &error);
            //cout << "# vegas warmup " << result_x << " +/- " << error_x << endl;
            int iter = 0;
            do
            {
                gsl_monte_vegas_integrate(&F, lower, upper, F.dim, MCINTPOINTS, global_rng, s, &result, &error);
                //cout << "# Vegas interation (" << complist[i] << ") " << result << " +/- " << error << " chisqr " << gsl_monte_vegas_chisq(s) << endl;
                iter++;
                
            } while (iter < 2 or ((fabs( gsl_monte_vegas_chisq(s) - 1.0) > 0.3 or std::abs(error/result) > VEGAS_RESULT_ACCURACY_TARGET) and iter < MAXITER_VEGAS));
            gsl_monte_vegas_free(s);
            
            if (std::abs(error/result) > VEGAS_RESULT_ACCURACY_TARGET and std::abs(result)>0)
            {
                cerr << "WARNING: Relative uncertainty (" << complist[i] << ") " << std::abs(error/result) << " at B=" << B <<", theta_B=" << theta_B << endl;
            }/*
            else {
                cerr << "OK: Relative uncertainty " << std::abs(error_x/result_x) << " at theta_B=" << theta_B << endl;
            }*/
        }
        
        results.push_back(result);
        
        
        
    }
   
    
    
    //if (std::abs(error/result) > MCINTACCURACY)
    //    cerr << "#MC integral failed, result " << result << " error " << error << endl;
    
    delete lower;
    delete upper;
    
    return results;
    
}



double Inthelperf_amplitude_z(double z, void* p);

double Is_point_charge(double B, void* p)
{
    Inthelper_amplitude *par = (Inthelper_amplitude*)p;
    double Z, A,  sqrts, mA;
    const double sqrt_aem = std::sqrt(1.0/137.0);
    const double MV  = 3.097; // JPsi
    if (KINEMATICS == LHC)
    {
        Z=82;
        sqrts=5020;
        A = 208;
        mA = 207.9766;
    }
    else if (KINEMATICS == RHIC)
    {
        Z = 79;
        A = 197;
        sqrts = 200;
        mA = 196.96656;
    }
    
    double y = -std::log(sqrts * par->xpom / MV);
    
    const double omega = MV/2.*std::exp(y);
    const double gamma = sqrts / (2.0*mA/A);
    
    return Z * sqrt_aem / M_PI * (omega/gamma) * gsl_sf_bessel_K1(B * omega/gamma);
    
    
}

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
    
    result = blen * std::exp(-imag*(q*b)) * subamp;
    
    double b_minus_Bv_len = (b-Bv).Len();
    //double b_plus_Bv_len = 0; // not used, save one sqrt() (b+Bv).Len();
    
    //if (b_plus_Bv_len < 1e-15 or b_minus_Bv_len < 1e-15)
    if(b_minus_Bv_len < 1e-15)
        return 0;
    
    double Is_b_minus_Bv;
    double Is_b_plus_Bv=0; // not used...
    
    if (NUCLEAR_FF == APPROXIMATIVE)
    {
        Is_b_minus_Bv =par->diffraction->GetIsInterpolator()->Evaluate(b_minus_Bv_len);
        //Is_b_plus_Bv = par->diffraction->GetIsInterpolator()->Evaluate(b_plus_Bv_len);
    }
    else if (NUCLEAR_FF == POINT_CHARGE)
    {
        Is_b_minus_Bv = Is_point_charge(b_minus_Bv_len, par);
       // Is_b_plus_Bv = Is_point_charge(b_plus_Bv_len);
    }
    else
    {
        cerr << "Unknown nuclear FF!" << endl;
        exit(1);
    }
    
    if (par->component == X)
    {
        cerr << "Support dropped for X comp" << endl; exit(1);
        //result *= (blen*cos_th_b - B*cos_th_B)/( b_minus_Bv_len ) * Is_b_minus_Bv
        //+ (blen*cos_th_b + B*cos_th_B)/( b_plus_Bv_len ) *  Is_b_plus_Bv * std::exp(-imag*(q*Bv));
    }
    else if (par->component == Y)
    {
        cerr << "Support dropped for Y comp" << endl; exit(1);
        //result *= (blen*sin_th_b - B*sin_th_B)/( b_minus_Bv_len ) * Is_b_minus_Bv
        //+ (blen*sin_th_b + B*sin_th_B)/( b_plus_Bv_len ) * Is_b_plus_Bv * std::exp(-imag*(q*Bv));
    }
    // In Farid's note there is (-iA) which is my subamp
    else if (par->component == M0)
    {
        result *= Is_b_minus_Bv / b_minus_Bv_len;
    }
    else if (par->component == M1x)
    {
        result *= cos_th_b * blen * Is_b_minus_Bv / b_minus_Bv_len;
    }
    else if (par->component == M1y)
    {
        result *= sin_th_b * blen * Is_b_minus_Bv / b_minus_Bv_len;
    }
    else
    {
        std::cerr << "Unknown component " << par->component << endl;
        exit(1);
    }
    double res;
    if (par->real_part)
        res = result.real();
    else
        res = result.imag();
    
    if (isnan(res))
        cerr << "What, result is NaN, b=" << b <<", r=" << r <<", z=" << z << endl;
    
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
    
    cerr << "Not used!" << endl; exit(1); 
    return 0;
    

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




double Diffraction::ScatteringAmplitudeRotationalSymmetry(double xpom, double Qsqr, double t, Polarization pol)
{
    return 0; // not implemented
}

double Diffraction::ScatteringAmplitudeRotationalSymmetryIntegrand(double xpom, double Qsqr, double t, double r, double b, double z, Polarization pol)
{
    return 0;
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
