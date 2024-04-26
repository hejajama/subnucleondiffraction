#ifndef VirtualPhoton_H
#define VirtualPhoton_H


#include "wave_function.hpp"
#include "qcd.hpp"
#include <iostream>
#include <string>
#include <vector>

/**
 * \class VirtualPhoton virtual_photon.hpp amplitudelib/virtual_photon.hpp
 * Virtual photon-q bar q overlap
 * 
 * Ref: Kowalski, Motyka and Watt, see arXiv: hep-ph/0606272v2
 * 
 * By default uses 3 light quarks (u,d,s)
 */
class VirtualPhoton : public WaveFunction {
    public:
        VirtualPhoton();
        
        /**
         * Overlap between q bar
         * q and transverse photon
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        double PsiSqr_T(double Qsqr, double r, double z);

        /**
         * Overlap between q bar q and longitudinal photon
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        double PsiSqr_L(double Qsqr, double r, double z);
        
        // Overlap wave functions integrated over z=[0,1]
        /**
         * Overlap between q bar q and transverse photon integrated over z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        double PsiSqr_T_intz(double Qsqr, double r);

        /**
         * Overlap between q bar q and longitudinal photon integrated over z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        double PsiSqr_L_intz(double Qsqr, double r);
        
        
        
        std::string GetParamString();
        
        /**
         * Select quar flavor
         *
         * If not set, sums over light quarks (u,d,s) with default mass 0.14 GeV.
         *
         * @param p quark flavor
         * @param mass quark mass, optional, default values used if not given
         */
        void SetQuark(Parton p, double mass=-1);

        std::string WaveFunctionType() { return "Photon"; }
        
        
    private:
        // Parameters
        std::vector<double> e_f;   // Quark charges
        std::vector<double> m_f;   // Quark masses (GeV)

        double Epsilon(double Qsqr, double z, int f);

        QCD qcd;
        
        
};

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic);



#endif
