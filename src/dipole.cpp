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
    gsl_root_fsolver_set(s, &f, 1e-6, 4);
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