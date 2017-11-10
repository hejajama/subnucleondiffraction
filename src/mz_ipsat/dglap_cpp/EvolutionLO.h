//==============================================================================
//  EvolutionLO.h
//
//  LO DGLAP evolution
//
//  Constructor initializes the Mellin momenta
//  and takes a valid instance of AlphaStrong.
//  Note: EvolutionLO takes ownership of
//        AlphaStrong.
//
//  The evolution functions xG() and alphasxG()
//  take as arguments:
//
//   x :       Bjorken momentum fraction
//   Q2:       Virtuality
//   mu0:      Input mu_r in GeV
//   coupling: 0 for decoupled evolution, 1 for coupled evolution
//   Ag:
//   lambdag:
//   As:
//   lambdas:
//
//  xG()       returns x*G(x, Q^2)
//  alphasxG() returns alpha_s(Q^2)*x*G(x,Q^2)
//
//==============================================================================
#ifndef EvolutionLO_h
#define EvolutionLO_h
#include <complex>
#include <iostream>
#include "AlphaStrong.h"

using namespace std;

class EvolutionLO {
public:
    EvolutionLO(AlphaStrong*);
    ~EvolutionLO();
    
    double xG(double x, double Q2, double mu0, int coupling, double Ag,
              double lambdag, double As, double lambdas);
    
    double alphasxG(double x, double Q2, double mu0, int coupling, double Ag,
                    double lambdag, double As, double lambdas);

private:
    void anom();
    void anCalc(complex<double>& qqi, complex<double>& qgf,
                  complex<double>& gqi, complex<double>& ggi,
                  complex<double>& ggf, complex<double>& xn);
    complex<double> psiFunction(complex<double>);
    complex<double> lngam(complex<double> X);
    complex<double> beta(complex<double>, complex<double>);
    void reno(complex<double>* fn, double alpq, int nmax, int coupling,
              double ag, double lambdag, double as, double lambdas);

private:
    double   mWN[137];
    complex<double>   mN[137];
    double   mC;
    complex<double> mCC;
    complex<double> mANS[137][6];
    complex<double> mAM[137][6];
    complex<double> mAP[137][6];
    complex<double> mAL[137][6];
    complex<double> mBE[137][6];
    complex<double> mAB[137][6];
    complex<double> mAC[137][6];
    double mALPS;
    double mALPC;
    double mALPB;
    double mALPT;
    double mFR2;
    double mMUR;
    double mASMUR;
    double mMC;
    double mMB;
    double mMT;
    
    AlphaStrong *mAlphaStrong;
};
#endif

