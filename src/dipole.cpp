/*
 * Dipole cross section with explicit impact parameter and r dependence
 * Abstract class
 * All nuclear (target) dependence is described by this class]
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "dipole.hpp"
#include "vector.hpp"
#include <string>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <cmath>

double DipoleAmplitude::Amplitude(double xpom, double q1[2], double q2[2] )
{
    return 0; // Not defined, shouldnt be here
}

void DipoleAmplitude::InitializeTarget()
{
    // Nothing to do if not overloaded
}

std::string DipoleAmplitude::InfoStr()
{
    return "Dipole modle InfoStr not implemented";
}

double DipoleAmplitude::Amplitude(double xpom, Vec q1, Vec q2)
{
    double quark[2] = {q1.GetX(), q1.GetY() };
    double antiquark[2] = {q2.GetX(), q2.GetY() };
    return Amplitude(xpom, quark, antiquark);
}

/*
 * Solve saturation scale
 * see dipole.hpp for definition
 *
 * Returns -1 if no satscale can be found (e.g. too large b in ipsat)
 */
struct SatscaleHelperDipole
{
    double xbj;
    double x,y;
    DipoleAmplitude* proton;
};
double SatscaleHelperfDipole(double r, void* p);
double DipoleAmplitude::SaturationScale(double xpom, Vec b)
{
    const int MAXITER = 100;
    const double ACC = 0.0001;
    gsl_function f;
    f.function = &SatscaleHelperfDipole;
    SatscaleHelperDipole par;
    par.xbj = xpom;
    par.x = b.GetX();
    par.y=b.GetY();
    par.proton = this;
    f.params = &par;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    
    // First find smallest r such that Q_s is between 0 and r
    // Here need to find smallest, as due to fluctuations there might be multiple solutions
    // to the equation
    double step = 0.1;
    double maxr = step;
    bool error=true;
    const double LIMITR = 100;
    for (maxr=step; maxr < LIMITR; maxr+=step)
    {
        if (SatscaleHelperfDipole(maxr, &par) > 0)
        {
            error=false;
            break;
        }
    }
    if (error)
    {
        //std::cerr << "Cant find saturation scale at b=" << b << ", SatscaleHelperfDipole(" << LIMITR <<")=" <<SatscaleHelperfDipole(LIMITR, &par) << std::endl;
        return -1;
    }
    
    gsl_root_fsolver_set(s, &f, 0, maxr);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ACC);
    } while (iter <= MAXITER and status == GSL_CONTINUE);
    
    double satr = gsl_root_fsolver_root(s);
    
    double qs = std::sqrt( 2.0) / satr;
    
    gsl_root_fsolver_free(s);
    return qs;
    
}


/*
 * Solve saturation scale defined as N(r^2=2/Qs^2) = 1 - exp(-1/2)
 */


double SatscaleHelperfDipole(double r, void* p)
{
    SatscaleHelperDipole *helper = (SatscaleHelperDipole*) p;
    Vec q1(helper->x-0.5*r, helper->y);
    Vec q2(helper->x+0.5*r, helper->y);
    // Now |q1-q2|=r, 0.5(q1+q2)=b
    double res =helper->proton->Amplitude(helper->xbj, q1, q2) - (1.0 - std::exp(-0.5));
    return res;
    

}


/*
 * If imaginary part of the amplitude is not implemented (e.g. ipsat), we set it to 0
 */
 double DipoleAmplitude::AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] )
{
    return 0;
}

double DipoleAmplitude::AmplitudeImaginaryPart(double xpom, Vec q1, Vec q2)
{
    double qq1[2] = {q1.GetX(), q1.GetY()};
    double qq2[2] = {q2.GetX(), q2.GetY()};
    return AmplitudeImaginaryPart(xpom, qq1, qq2);
}

void DipoleAmplitude::SetSkewedness(bool s)
{
    std::cerr << "DipoleAmplitude does not have SetSkewedness implemented, ignoring..." << std::endl;
}

double DipoleAmplitude::Density(Vec b)
{
    std::cerr << "DipoleAmplitude does not have Density implemented, ignoring..." << std::endl;
    return 0;
}
