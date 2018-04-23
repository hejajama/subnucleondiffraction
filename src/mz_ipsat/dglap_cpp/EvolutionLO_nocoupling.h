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
#ifndef EvolutionLO_gluon_h
#define EvolutionLO_gluon_h
#include <complex>
#include <iostream>
#include "AlphaStrong.h"

using namespace std;

class EvolutionLO_gluon {
public:
    EvolutionLO_gluon(AlphaStrong*);
    ~EvolutionLO_gluon();
    
    double xG(double x, double Q2, double mu0, int coupling, double Ag,
              double lambdag, double As, double lambdas);
    
    double alphasxG(double x, double Q2, double mu0, int coupling, double Ag,
                    double lambdag, double As, double lambdas);
    
    void generateLookupTable( double mu0, int coupling, double Ag,
                             double lambdag, double As, double lambdas, int nx = 200, int nq2 = 200);
    void useLookupTable(bool) ;
    bool lookupTableIsUsed() const;
    double xG_Interpolator(double x, double Q2);
    double luovi( double f[4], double arg[4], double z);

private:
    void anom();
    void anCalc( complex<double>& ggi,
                  complex<double>& ggf, complex<double>& xn);
    complex<double> psiFunction(complex<double>);
    complex<double> lngam(complex<double> X);
    complex<double> beta(complex<double>, complex<double>);
    void reno(complex<double>* fn, double alpq, int nmax, int coupling,
              double ag, double lambdag, double as, double lambdas);

private:
    static EvolutionLO_gluon* mInstance;
    double   mWN[137];
    complex<double>   mN[137];
    double   mC;
    complex<double> mCC;
    complex<double> mAP[137][6];
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
    
    bool mLookupTableIsFilled;
    bool mUseLookupTable;
    unsigned int mNumberOfNodesInX;
    unsigned int mNumberOfNodesInQ2;
    double mTableMinX;
    double mTableMaxX;
    double mTableMinQ2;
    double mTableMaxQ2;
    double **mLookupTable;
    
    AlphaStrong *mAlphaStrong;
};

inline void EvolutionLO_gluon::useLookupTable(bool val) {mUseLookupTable = val;}

inline bool EvolutionLO_gluon::lookupTableIsUsed() const
{
    return mUseLookupTable && mLookupTableIsFilled;
}



#endif

