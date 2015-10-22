/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "wave_function.hpp"
#include <iostream>

WaveFunction::WaveFunction()
{
    mode=VM_MODE_TOT;
}

/*
 * PsiSqr_intz
 * Returns |\Psi_L|^2, |\Psi_T|^2 or |\Psi_L|^2 + |\Psi_T|^2 depending
 * on the specified mode
 */
REAL WaveFunction::PsiSqr_intz(REAL Qsqr, REAL r)
{
    switch (mode) 
    {
        case VM_MODE_TOT:
            std::cerr << "Can't calculate total wave function overlap in "
             << "WaveFunction::PsiSqr_intz!" << std::endl;
            break;
        case VM_MODE_L:
            return PsiSqr_L_intz(Qsqr, r);
            break;
        case VM_MODE_T:
            return PsiSqr_T_intz(Qsqr, r);
            break;
        default:
            std::cerr << "WaveFunctioN::PsiSqr_intz: Unknown mode " << mode 
                << std::endl;
     }
    return 0; // Shouldnt be here!
}  

REAL WaveFunction::PsiSqr_tot(REAL Qsqr, REAL r, REAL z)
{
    return PsiSqr_T(Qsqr,r,z)+PsiSqr_L(Qsqr,r,z);
}

REAL WaveFunction::PsiSqr_tot_intz(REAL Qsqr, REAL r)
{
    return PsiSqr_T_intz(Qsqr,r)+PsiSqr_L_intz(Qsqr,r);
}

void WaveFunction::SetMode(int m)
{
    mode=m;
}



// eps(z,Q,r) = sqrt(z(1-z)*Q^2 + m^2)
double epsfunsqr(double z, double Qsqr, double msqr)
{
    return z*(1.0-z)*Qsqr + msqr;
}

double epsfun(double z, double Qsqr, double msqr) {
    return std::sqrt(epsfunsqr(z,Qsqr,msqr));
}
