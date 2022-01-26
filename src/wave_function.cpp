/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "wave_function.hpp"
#include <iostream>

WaveFunction::WaveFunction()
{

}


double WaveFunction::PsiSqr_tot(double Qsqr, double r, double z)
{
    return PsiSqr_T(Qsqr,r,z)+PsiSqr_L(Qsqr,r,z);
}

double WaveFunction::PsiSqr_tot_intz(double Qsqr, double r)
{
    return PsiSqr_T_intz(Qsqr,r)+PsiSqr_L_intz(Qsqr,r);
}

double WaveFunction::MesonMass()
{
    return 0;
}


std::string WaveFunction::WaveFunctionType()
{
	return "Not specified";
}

// eps(z,Q,r) = sqrt(z(1-z)*Q^2 + m^2)
double epsfunsqr(double z, double Qsqr, double msqr)
{
    return z*(1.0-z)*Qsqr + msqr;
}

double epsfun(double z, double Qsqr, double msqr) {
    return std::sqrt(epsfunsqr(z,Qsqr,msqr));
}


