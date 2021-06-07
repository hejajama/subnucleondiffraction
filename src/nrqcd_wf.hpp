#ifndef NRQCDWF_H
#define NRQCDWF_H

/*
 * Overlap between the photon and the vector meson wave functions
 * NRQCD parametrization
 * T. Lappi, H. Mäntysaari, J. Penttala, arXiv:2006.02830
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2021
 */

#include <tools/tools.hpp>
#include <amplitudelib/wave_function.hpp>
#include <iostream>
#include <string>


class NRQCD_WF : public WaveFunction {
    public:
    NRQCD_WF(double A_=0.213, double B_=-0.0157, double ef_=2.0/3.0, double mc_=1.4, double MV_=3.097);
        
        
        // Overlap wave functions
        double PsiSqr_T(double Qsqr, double r, double z);
        double PsiSqr_L(double Qsqr, double r, double z);
        
        // Overlap wave functions integrated over z=[0,1]
        // NB! Normalization factor 1/(4*Pi) is included here!
        // So calculates \int dz/(4*Pi) PsiSqr_T/L exp(i (z-1/2) r.Delta)
        // The result depends on angle between r,Delta!
        double PsiSqr_T_intz(double Qsqr, double r, double Delta, double phi_r_Delta);
        double PsiSqr_L_intz(double Qsqr, double r, double Delta, double phi_r_Delta);
    
        double PsiSqr_T_intz(double Qsqr, double r) { return 0; }
        double PsiSqr_L_intz(double Qsqr, double r) { return 0; }
        
        std::string GetParamString();

        double MesonMass();
    
        std::string WaveFunctionType() { return "NRQCD"; }

    
        
    private:
        // Parameters
    double A,B;
    double MV,mc;
    double ef;
};

std::ostream& operator<<(std::ostream& os, NRQCD_WF& ic);

#endif // NRQCDWF_H
