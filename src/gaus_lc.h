#ifndef GausLC_H
#define GausLC_H

/*
 * Overlap between the photon and the vector meson wave functions
 * Boosted gaussian version
 * Ref: Kowalski, Motyka and Watt, see arXiv: hep-ph/0606272v2
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include <tools/tools.hpp>
#include <amplitudelib/wave_function.hpp>
#include <iostream>
#include <string>

typedef double REAL;

class GausLC : public WaveFunction {
    public:
        GausLC(REAL e_f_, REAL N_T_, REAL N_L_, REAL R_T_, REAL R_L_, 
                REAL m_f_, REAL M_V_, int delta_);
        GausLC(std::string file);
        
        // Overlap wave functions
        REAL PsiSqr_T(REAL Qsqr, REAL r, REAL z);
        REAL PsiSqr_L(REAL Qsqr, REAL r, REAL z);
        
        // Overlap wave functions integrated over z=[0,1]
        // NB! Normalization factor 1/(4*Pi) is included here!
        // So calculates \int dz/(4*Pi) PsiSqr_T/L
        REAL PsiSqr_T_intz(REAL Qsqr, REAL r);
        REAL PsiSqr_L_intz(REAL Qsqr, REAL r);
        
        // Some functions required by PsiSqr_T/L
        REAL Psi_T(REAL r, REAL z);
        REAL Psi_L(REAL r, REAL z);
        REAL Psi_T_DR(REAL r, REAL z); // \partial_r Psi_T
        REAL Psi_L_DR(REAL r, REAL z); 
        REAL Psi_L_D2R(REAL r, REAL z); // \partial^2_r Psi_L

        REAL MesonMass();
        
        std::string GetParamString();


        double QuarkMass() { return m_f; }
        int GetDelta() { return delta; } 
        
        
    private:
        // Parameters
        REAL e_f;   // Charge
        REAL N_T, N_L; // Constants for Psi_T and Psi_L
        REAL R_T, R_L;
        REAL m_f,M_V;
        int delta;  // There are two different longitudial wave functions,
                    // delta=0 or delta=1
        
        
};

std::ostream& operator<<(std::ostream& os, GausLC& ic);

#endif // GausLC_H
