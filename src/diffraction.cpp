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
#include "gauss_boost.hpp"

using namespace std;

#include <complex>


Diffraction::Diffraction(DipoleAmplitude& dipole_, WaveFunction& wavef_)
{
    dipole=&dipole_;
    wavef=&wavef_;
    num_of_averages = 1;
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
    double res = std::log(par->diff->ScatteringAmplitudeRotationalSymmetry(x, par->Qsqr, par->t, par->pol));
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
    DVCS_COMPONENT comp;
};

double Inthelperf_amplitude_mc( double *vec, size_t dim, void* par);

double Diffraction::ScatteringAmplitude(double xpom, double Qsqr, double t, DVCS_COMPONENT comp)
{
    Inthelper_amplitude helper;
    helper.diffraction = this;
    helper.xpom = xpom;
    helper.Qsqr = Qsqr;
    helper.t = t;
    helper.comp=comp;

    
    
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
    upper[0] = 15*5.068 ; // Max b
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
    
    if (!FACTORIZE_ZINT)
        z = vec[4];
        
    return helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z, helper->comp);
    
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
    
    return 0; //helper->diffraction->ScatteringAmplitudeIntegrand(helper->xpom, helper->Qsqr, helper->t, helper->r, helper->theta_r, helper->b, helper->theta_b, z);
}

double Diffraction::ScatteringAmplitudeIntegrand(double xpom, double Qsqr, double delta, double r, double theta_r, double b, double theta_b, double z, DVCS_COMPONENT comp)
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
    if (FACTORIZE_ZINT)
    {
        cerr << "Facrtorized zint not supported!" << endl;
        exit(1);
    }
    
    double qx = bx + z*rx; double qy = by + z*ry;
    double qbarx = bx - (1.0-z)*rx; double qbary = by - (1.0-z)*ry;
    
    z = tmpz;

    std::complex<double> res = 0;
    
    res = r * b;  // r and b from Jacobians
    // Compared to subnulceondiff code, no factor 2
   
    // For VM, need scalar part. Assume BG here
    BoostedGauss *BG = (BoostedGauss*)wavef;
    double MV = BG->MesonMass();
 
    //double delta = std::sqrt(t);
    
    double x1[2] = {qx,qy};
    double x2[2] = {qbarx, qbary};
    double amp_real = dipole->Amplitude(xpom, x1, x2 );
    //double amp_imag = 0; // todo: check imaginary part, should vanihs... dTipole->AmplitudeImaginaryPart(xpom, x1, x2);
    double amp_imag = 0; //dipole->AmplitudeImaginaryPart(xpom, x1, x2);
    std::complex<double> amp(amp_real, amp_imag);
    //amp = amp.real();   // Disable possible imag part for now
    //
    // DEGBUG TEST
    //double cb=0.0;
    //amp = 1.0 - std::exp(-r*r*1*1/4.0*std::exp(-b*b/8.0) * (1.0 - cb*(0.5 - SQR(std::cos(theta_r - theta_b)) ) ) );
    
    double phasedelta = (2.0*z - 1.0)/2.0*delta*r*std::cos(theta_r);
        
    std::complex<double> imag(0,1);
    double mf = 0.14;
    double vm_phi_l_wf=0; // [Mv + (mf^2 - Nabla^2)/(Mv*z(-1z))] * phi_L
    if (comp == VM_LL or comp==VM_TT or comp==VM_TTflip or comp==VM_LT or comp==VM_TL)
    {
        mf = BG->QuarkMass();
        vm_phi_l_wf = (MV + mf*mf/(MV*z*(1.-z))) * BG->Psi_L(r,z) - 1.0/(MV*z*(1.-z)) * ( 1.0/r*BG->Psi_L_DR(r,z) + BG->Psi_L_D2R(r,z) );
    }

    double Q = std::sqrt(Qsqr);
    // Eqs from Farids note (60)-(64)
    // Without prefactor (4 Nc qf^2)^2/(2pi)^2

    double epscale = std::sqrt(z*(1.0-z)*Qsqr + mf*mf); 
    double epscale2 = mf; // Q'^2=0


    // NOTE: Compared to our paper, an overall minus sign is missing. But this sign has no effect on cross sections
   if (comp == TT)
       res *= std::exp(-imag*( b*delta*std::cos(theta_b) + phasedelta)) * (
             (z*z + SQR(1.0-z)) * epscale * gsl_sf_bessel_K1(epscale*r) * epscale2*gsl_sf_bessel_K1(epscale2*r)
             + mf*mf*gsl_sf_bessel_K0(epscale*r)*gsl_sf_bessel_K0(epscale2*r) )
        * amp;
   else if (comp == LL)
        res = 0;
    else if (comp == TTflip)
        res*= -2.*std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta - 2.0*theta_r)) * z*(1.0-z)* epscale * gsl_sf_bessel_K1(epscale*r) * epscale2*gsl_sf_bessel_K1(epscale2*r) * amp;
   else if (comp == LT)
        res *= -std::sqrt(2) * imag * std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta - theta_r)) * z*(1.0-z)*(2.0*z-1.)*std::sqrt(Qsqr) * gsl_sf_bessel_K0(epscale*r)* epscale2*gsl_sf_bessel_K1(epscale2*r) * amp;
    else if (comp == TL)
        res = 0;
    else if (comp == VM_LL)
    {
        res *= 2.*std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta)) * SQR(z*(1.-z))/(z*(1.-z)) * Q * gsl_sf_bessel_K0(epscale*r) * vm_phi_l_wf * amp;    
    }
    else if (comp == VM_TT) 
    {
        res *= - std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta)) * ( 
                (SQR(z)+SQR(1.-z))/(z*(1.-z)) * epscale * gsl_sf_bessel_K1(epscale*r) * BG->Psi_T_DR(r,z) * amp
                - 1./(z*(1.-z)) * mf*mf*gsl_sf_bessel_K0(epscale*r)*BG->Psi_T(r,z) * amp );
    }
    else if (comp == VM_TTflip)
    {
        res *= 2.*std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta - 2.0*theta_r)) * epscale * gsl_sf_bessel_K1(epscale*r) * BG->Psi_T_DR(r,z)*amp;
    }
    else if (comp == VM_LT)
    {
        res *= imag*std::sqrt(2) * std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta - theta_r)) * (2.0*z-1.0) * Q * gsl_sf_bessel_K0(epscale*r) * BG->Psi_T_DR(r,z)*amp;
   }
    else if (comp == VM_TL)
    {
        // Sign fixed 
        res *= imag/std::sqrt(2) * std::exp(-imag*(b*delta*std::cos(theta_b) + phasedelta + theta_r)) * (2.0*z-1.) * epscale * gsl_sf_bessel_K1(epscale*r) * vm_phi_l_wf * amp;
    }
    else
        {
        cerr << "Unknown component!" << endl; exit(1);
    }
    // dz measure
    res /= 4.0*M_PI;
    ///NOTE: in main.cpp I multiply by 4pi!

    if (REAL_PART)
        return res.real();
    else
        return res.imag();
    

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
    
    return 0;;
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
    return 0;
    //return par->diffraction->ScatteringAmplitudeRotationalSymmetryIntegrand(par->xpom, par->Qsqr, par->t, par->r, par->b, z, par->polarization);
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
