#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
/** @class */
/*
 * AmplitudeLib
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2014
 */

#include <string>
#include <cmath>
/**
 * Virtual class from which different wave functions are inherited
 *
 * @see VirtualPhoton
 */
class WaveFunction{
    public:
        WaveFunction();
        virtual ~WaveFunction() = default;
        
        /**
         * Transverse overlap
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        virtual double PsiSqr_T(double Qsqr, double r, double z) = 0;

        /**
         * Longitudinal overlap
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        virtual double PsiSqr_L(double Qsqr, double r, double z) = 0;
        /**
         * Transverse overlap integrated over mom. fraction z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        virtual double PsiSqr_T_intz(double Qsqr, double r) = 0;
        /**
         * Longitudinal overlap integrated over mom. fraction z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        virtual double PsiSqr_L_intz(double Qsqr, double r) = 0;
        /**
         * Return string containing info on parameters
         */
        virtual std::string GetParamString()=0;
        /**
         * Longitudinal overlap^2 + transverse overlap^2
         */
        double PsiSqr_tot(double Qsqr, double r, double z);
        /**
         * Longitudinal overlap^2 + transverse overlap^2
         * integrated over z.
         */
        double PsiSqr_tot_intz(double Qsqr, double r);
    
        /*
         * Vector meson mass
         */
        virtual double MesonMass();

		virtual std::string WaveFunctionType();
};

// Helpers
double epsfunsqr(double z, double Qsqr, double msqr);
double epsfun(double z, double Qsqr, double msqr);
#endif  // WAVE_FUNCTION_H
