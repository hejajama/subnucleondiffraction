//--------------------------------------------------------------------
//
//   Class AlphaStrong
//
//   This is a rewrite of the original alphaS.f code originally
//   written in F77. The documentation is largely left as provided
//   in the Fortran version but includes the updated variable and
//   function names. Data originally stored in common
//   blocks and several static variables turned into data members.
//   Some relics of the original F77 code remained such as array
//   indexing starting at 1 (arrays are size +1 here).
//
//   Original author: Graeme Watt
//   C++ version: Thomas Ullrich (BNL)
//
//--------------------------------------------------------------------
//   Stand-alone code for alpha_s cannibalised (with permission)
//   from Andreas Vogt's QCD-PEGASUS package (hep-ph/0408244).
//   The running coupling alpha_s is obtained at N^mLO (m = 0,1,2,3)
//   by solving the renormalisation group equation in the MSbar scheme
//   by a fourth-order Runge-Kutta integration.  Transitions from
//   n_f to n_f+1 flavours are made when the factorisation scale
//   mu_f equals the pole masses m_h (h = c,b,t).  At exactly
//   the thresholds m_{c,b,t}, the number of flavours n_f = {3,4,5}.
//   The top quark mass should be set to be very large to evolve with
//   a maximum of five flavours.  The factorisation scale mu_f may be
//   a constant multiple of the renormalisation scale mu_r.  The input
//   factorisation scale mu_(f,0) should be less than or equal to
//   the charm quark mass.  However, if it is greater than the
//   charm quark mass, the value of alpha_s at mu_(f,0) = 1 GeV will
//   be found using a root-finding algorithm.
//
//   Example of usage.
//   The constructor takes the following arguments:
//
//    order = 2                 // perturbative order (N^mLO,m=0,1,2,3)
//    fr2 = 1.                  // ratio of mu_f^2 to mu_r^2
//    mur = 1.                  // input mu_r in GeV
//    asmur = 0.5               // input value of alpha_s at mu_r
//    mc = 1.4                  // charm quark mass
//    mb = 4.75                 // bottom quark mass
//    mt = 1.e10                // top quark mass
//
//    AlphaStrong myAlphaS(order, fr2, mur, asmur, mc, mb, mt);
//
//   Then get alpha_s at a renormalisation scale mu_r with:
//
//    mu_r = 100.;              // renormalisation scale in GeV
//    double result = myAlphaS.value(mu_r);
//
//--------------------------------------------------------------------
#ifndef AlphaStrong_h
#define AlphaStrong_h
#include <cmath>

using namespace std;

class AlphaStrong {
public:
    AlphaStrong(int order, double fr2, double mur, double asmur, double mc, double mb, double mt);

    double value(double);
    
    int order() const;            // used order
    double massCharm() const;     // used charm quark mass
    double massBottom() const;    // used bottom quark mass
    double massTop() const;       // used bottom quark mass
    double ratioFR2() const;      // used ratio of mu_f^2 to mu_r^2
    double alphasAtMuR() const;   // input value of alpha_s at mu_r

private:
    double findR0(double asi);
    double asnf1 (double asnf, double logrh, int nf);
    double as(double r2, double r20, double as0, int nf);
    double funBeta1(double, int);
    double funBeta2(double, int);
    double funBeta3(double, int);
    double rootFinder(double a0, double b0, double eps, int maxf, int mode);
    double sign(double, double);

    void evolution (double mc2, double mb2, double mt2);
    void initR0(int order, double fr2, double r0, double asi, double mc, double mb, double mt);
    void betaFunction();

private:
    double mFR2;
    double mScale;
    double mAsScale;
    double mMassC;
    double mMassB;
    double mMassT;
    double mR0;
    double mOrder;
    double mZeta[7];
    double mCF;
    double mCA;
    double mTR;
    double mAS0;
    double mM20;
    double mLogFR;
    double mASC;
    double mM2C;
    double mASB;
    double mM2B;
    double mAST;
    double mM2T;
    double mBeta0[7];
    double mBeta1[7];
    double mBeta2[7];
    double mBeta3[7];
    double mCCMCoefficients[4][4];
    double mCCMCoefficientsI30;
    double mCCMCoefficientsI31;
    double mCCMCoefficientsF30;
    double mCCMCoefficientsF31;
    int mWasCalled;
    int mPertubativeOrder;
    int mIntegrationSteps;
    int mVarFlavourNumScheme;
    int mNumFlavorsFFNS;
};

inline int AlphaStrong::order() const {return mOrder;}

inline double AlphaStrong::massCharm() const  {return mMassC;}

inline double AlphaStrong::massBottom() const {return mMassB;}

inline double AlphaStrong::massTop() const {return mMassT;}

inline double AlphaStrong::ratioFR2() const {return mFR2;}

inline double AlphaStrong::alphasAtMuR() const {return mAsScale;}

inline double AlphaStrong::sign(double a, double b) {
    if (b < 0)
        return -fabs(a);
    else
        return fabs(a);
}
#endif
