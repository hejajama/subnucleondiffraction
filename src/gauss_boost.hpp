#ifndef BoostedGaus_H
#define BoostedGaus_H

/*
 * Overlap between the photon and the vector meson wave functions
 * Boosted Gaussian version, ref.
 * Kowalski, Motyka and Watt, see arXiv: hep-ph/0606272v2
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2014
 */

#include <tools/tools.hpp>
#include <amplitudelib/wave_function.hpp>
#include <iostream>
#include <string>

typedef double REAL;

/* Overlap functions are taken from an article written by Kowalski, Motyka 
 * and Watt, see arXiv: hep-ph/0606272v2
 */

class BoostedGauss : public WaveFunction {
    public:
        BoostedGauss(REAL e_f_, REAL N_T_, REAL N_L_, REAL R_T_, 
                REAL m_f_, REAL M_V_, int delta_);
        BoostedGauss(std::string file);
        
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
        
        std::string GetParamString();

        REAL MesonMass();
    
        void ZLimit(double zlim) { MINZ = zlim; MAXZ = 1.0-zlim; }
    
        
    private:
        // Parameters
        REAL e_f;
        REAL N_T, N_L; // Constants for Psi_T and Psi_L
        REAL R;
        REAL m_f,M_V;
        REAL alpha; // needed for higher exited states like 2s
        int S;      // S=1 is J/Psi etc.
        int delta;  // There are two different longitudial wave functions,
                    // delta=0 or delta=1
        double MINZ,MAXZ;
};

std::ostream& operator<<(std::ostream& os, BoostedGauss& ic);

#endif // BoostedGauss_H
