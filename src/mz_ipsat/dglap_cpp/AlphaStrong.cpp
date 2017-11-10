#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "AlphaStrong.h"

AlphaStrong::AlphaStrong(int order, double fr2, double mur, double asmur, double mc, double mb, double mt)
{
    //   order = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
    //   fr2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
    //   mur = input renormalisation scale (in GeV) for alpha_s.
    //   asmur = input value of alpha_s at the renormalisation scale mur.
    //   mc,mb,mt = heavy quark masses in GeV.
    
    const double eps = 1e-10;
    const int maxf = 10000;
    const int mode = 1;
    
    mWasCalled = 0;
    
    double r0, asi, a, b;
    
    mMassC = mc;
    mMassB = mb;
    mMassT = mt;
    mR0 = 1./sqrt(fr2);
    mOrder = order;
    mFR2 = fr2;
    mScale = mur;
    mAsScale = asmur;

    if (mur*sqrt(fr2) <= mc) {
        r0 = mur;
        asi = asmur;
    }
    else {                     // Solve for alpha_s at R0 = 1 GeV.
        a = 0.02;              // lower bound for alpha_s(R0)
        b = 2.00;              // upper bound for alpha_s(R0)
        r0 = mR0;
        //   Now get alpha_s(R0) corresponding to alpha_s(MUR).
        asi = rootFinder(a, b, eps, maxf, mode);
    }
    
    initR0(order, fr2, r0, asi, mc, mb, mt);
}

double AlphaStrong::findR0(double asi)
{
    //
    //  Function for dzerox
    //
    initR0(mOrder, mFR2, mR0, asi, mMassC, mMassB, mMassT);
    return value(mScale) - mAsScale; // solve equal to zero
}

void AlphaStrong::initR0(int order, double fr2, double R0, double ASI, double MC, double MB, double MT)
{
    //   order = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
    //   fr2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
    //   R0 = input renormalisation scale (in GeV) for alphas_s.
    //   ASI = input value of alpha_s at the renormalisation scale R0.
    //   MC,MB,MT = heavy quark masses in GeV.
    //   Must have R0*sqrt(fr2) <= MC to call this function.
    
    //
    //   QCD colour factors
    //
    mCA = 3.;
    mCF = 4./3.;
    mTR = 0.5;
    
    //
    //   The lowest integer values of the Zeta function
    //
    mZeta[1] = 0.57721566490153;
    mZeta[2] = 1.644934066848226;
    mZeta[3] = 1.202056903159594;
    mZeta[4] = 1.082323233711138;
    mZeta[5] = 1.036927755143370;
    mZeta[6] = 1.017343061984449;

    mVarFlavourNumScheme = 1;    // variable flavour-number scheme (VFNS)
    //mVarFlavourNumScheme = 0;  // fixed flavour-number scheme (FFNS)
    mNumFlavorsFFNS = 4;         // number of flavours for FFNS
    mPertubativeOrder = order;    // perturbative order of alpha_s
    mIntegrationSteps = 20;      // num. steps in Runge-Kutta integration
    double R20 = R0*R0;          // input renormalisation scale
    double MC2 = MC*MC;          // mu_f^2 for charm threshold
    double MB2 = MB*MB;          // mu_f^2 for bottom threshold
    double MT2 = MT*MT;          // mu_f^2 for top threshold
    mLogFR = log(fr2);           // log of ratio of mu_f^2 to mu_r^2
    mM20 = R20 * fr2;            // input factorisation scale

    //
    //   Stop some nonsense
    //
    if ( (mVarFlavourNumScheme == 0) && (mNumFlavorsFFNS < 3) ) {
        cout << "AlphaStrong::initR0(): Wrong flavour number for FFNS evolution. STOP." << endl;
        exit(1);
    }
    if ( (mVarFlavourNumScheme == 0) && (mNumFlavorsFFNS > 5) ) {
        cout << "AlphaStrong::initR0(): Wrong flavour number for FFNS evolution. STOP" << endl;
        exit(1);
    }
    
    if ( mPertubativeOrder > 3 ) {
        cout << "AlphaStrong::initR0(): Specified order in a_s too high ("
             << "mM20 = " << mM20 << ", MC2 = " << MC2 << "). STOP" << endl;
        exit(1);
    }
    
    if ( (mVarFlavourNumScheme != 0) && (fr2 > 4.001) ) {
        cout << "AlphaStrong::initR0(): Too low mu_r for VFNS evolution. STOP" << endl;
        exit(1);
    }
    
    if ( (mVarFlavourNumScheme == 1) && (mM20 > MC2) ) {
        cout << "AlphaStrong::initR0(): Too high mu_0 for VFNS evolution ("
             << "mM20 = " << mM20 << ", MC2 = " << MC2 << "). STOP" << endl;
        exit(1);
    }
    
    if ( (ASI > 2.) || (ASI < 2.e-2) ) {
        cout << "AlphaStrong::initR0(): alpha_s out of range. STOP" << endl;
        exit(1);
    }
    
    if ( (mVarFlavourNumScheme == 1) && (MC2 > MB2) ) {
        cout << "AlphaStrong::initR0(): Wrong charm-bottom mass hierarchy. STOP" << endl;
        exit(1);
    }
    if ( (mVarFlavourNumScheme == 1) && (MB2 > MT2) ) {
        cout << "AlphaStrong::initR0(): Wrong bottom-top mass hierarchy. STOP" << endl;
        exit(1);
    }
    
    betaFunction();
    
    // Keep a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
    mAS0 = ASI / (4*M_PI);
    
    if (mVarFlavourNumScheme != 0) {
        evolution (MC2, MB2, MT2);
    }
}

double AlphaStrong::value(double MUR)
{
    double R2 = MUR*MUR;
    double M2 = R2 * exp(mLogFR);
    int NF;
    double R20, ASI, ASF, R2T, R2B, R2C;
    
    if (mVarFlavourNumScheme == 0) {
        //
        //   Fixed number of flavours
        //
        NF  = mNumFlavorsFFNS;
        R20 = mM20 * R2/M2;
        ASI = mAS0;
        ASF = as(R2, R20, mAS0, NF);
    }
    else {
        //
        //   Variable number of flavours
        //
        if (M2 > mM2T) {
            NF = 6;
            R2T = mM2T * R2/M2;
            ASI = mAST;
            ASF = as(R2, R2T, mAST, NF);
        }
        else if (M2 > mM2B) {
            NF = 5;
            R2B = mM2B * R2/M2;
            ASI = mASB;
            ASF = as(R2, R2B, mASB, NF);
        }
        else if  (M2 > mM2C) {
            NF = 4;
            R2C = mM2C * R2/M2;
            ASI = mASC;
            ASF = as(R2, R2C, mASC, NF);
        }
        else {
            NF = 3;
            R20 = mM20 * R2/M2;
            ASI = mAS0;
            ASF = as(R2, R20, mAS0, NF);
        }
    }
    //
    //   Final value of alpha_s
    //
    double result = 4.*M_PI*ASF;
    return result;
}

double AlphaStrong::asnf1(double ASNF, double LOGRH, int NF)
{
    //
    //   The threshold matching of the QCD coupling in the MS(bar) scheme,
    //    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours
    //    up to order a_s^4 (NNNLO).
    //
    //   The value  ASNF  of a_s for NF flavours at the matching scale, the
    //    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
    //    of the heavy quark -- and  NF  are passed as arguments to the
    //    function  asnf1().  The order of the expansion  NAORD  (defined as
    //    the 'n' in N^nLO) is provided as data member(s).
    //
    //   The matching coefficients are inverted from Chetyrkin, Kniehl and
    //    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
    //    factors have been hard-wired in these results. The lowest integer
    //    values of the Zeta function are stored as data member.
    //

    //
    //   The coupling-constant matching coefficients (CMC's) up to NNNLO
    //   (calculated and saved in the first call of this routine)
    //
    if (mWasCalled != 1) {
        mCCMCoefficients[1][0] =  0.;
        mCCMCoefficients[1][1] =  2./3.;
        mCCMCoefficients[2][0] = 14./3.;
        mCCMCoefficients[2][1] = 38./3.;
        mCCMCoefficients[2][2] =  4./9.;
        mCCMCoefficientsI30 = + 80507./432. * mZeta[3] + 58933./1944. + 128./3. * mZeta[2] * (1.+ log(2.)/3.);
        mCCMCoefficientsF30 = - 64./9. * (mZeta[2] + 2479./3456.);
        mCCMCoefficientsI31 =   8941./27.;
        mCCMCoefficientsF31 = - 409./27.;
        mCCMCoefficients[3][2] = 511./9.;
        mCCMCoefficients[3][3] = 8./27.;
        mWasCalled = 1;
    }
    
    //
    //   The N_f dependent CMC's, and the alpha_s matching at order NAORD
    //
    mCCMCoefficients[3][0] = mCCMCoefficientsI30 + NF * mCCMCoefficientsF30;
    mCCMCoefficients[3][1] = mCCMCoefficientsI31 + NF * mCCMCoefficientsF31;
    
    double result = ASNF;
    if (mPertubativeOrder == 0) return result;
    double ASP = ASNF;
    
    for (int K1 = 1; K1 <= mPertubativeOrder; K1++) {
        ASP = ASP * ASNF;
        double LRHP = 1.;
        for (int K2 = 0; K2 <= K1; K2++) {
            result = result + ASP * mCCMCoefficients[K1][K2] * LRHP;
            LRHP = LRHP * LOGRH;
        }
    }
    return result;
}

void AlphaStrong::evolution (double MC2, double MB2, double MT2)
{
    //
    //   The function evolution() for the evolution of  a_s = alpha_s/(4 pi)
    //    from a three-flavour initial scale to the four- to six-flavour
    //    thresholds (identified with the squares of the corresponding quark
    //    masses).  The results are kept as data member.
    //
    //   The input scale  mM20 = mu_(f,0)^2  and the corresponding value
    //    mAS0  of a_s  are provided as data member as is the fixed scale
    //    logarithm mLogFR = ln (mu_f^2/mu_r^2).  The alpha_s
    //    matching is done by the function asnf1.
    //

    //
    //   Coupling constants at and evolution distances to/between thresholds
    //
    double R20 = mM20 * exp(-mLogFR);
    //
    //   Charm
    //
    mM2C  = MC2;
    double R2C  = mM2C * R20/mM20;
    double ASC3 = as(R2C, R20, mAS0, 3);
    // double SC   = log (mAS0 / ASC3);     // not used
    mASC  = asnf1 (ASC3, -mLogFR, 3);
    //
    //   Bottom
    //
    mM2B  = MB2;
    double R2B  = mM2B * R20/mM20;
    double ASB4 = as(R2B, R2C, mASC, 4);
    // double SB   = log (mASC / ASB4);     // not used
    mASB  = asnf1 (ASB4, -mLogFR, 4);
    //
    //   Top
    //
    mM2T  = MT2;
    double R2T  = mM2T * R20/mM20;
    double AST5 = as(R2T, R2B, mASB, 5);
    // double ST   = log (mASB / AST5);   // not used
    mAST  = asnf1 (AST5, -mLogFR, 5);
}

//
//   The beta functions funBeta'n' at N^nLO for n = 1, 2, and 3
//
double AlphaStrong::funBeta1(double A, int NF) {return - A*A * (mBeta0[NF] + A *   mBeta1[NF]);}

double AlphaStrong::funBeta2(double A, int NF) {return - A*A * (mBeta0[NF] + A * (mBeta1[NF] + A * mBeta2[NF]));}

double AlphaStrong::funBeta3(double A, int NF) {return - A*A * (mBeta0[NF] + A * (mBeta1[NF] + A * (mBeta2[NF] + A * mBeta3[NF])));}

double AlphaStrong::as(double R2, double R20, double AS0, int NF)
{
    //
    //   The running coupling of QCD,
    //
    //         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
    //
    //    obtained by integrating the evolution equation for a fixed number
    //    of massless flavours  NF.  Except at leading order (LO),  AS  is
    //    obtained using a fourth-order Runge-Kutta integration.
    //
    //   The initial and final scales  R20  and  R2,  the value  AS0  at
    //    R20, and  NF  are passed as function arguments.  The coefficients
    //    of the beta function up to  a_s^5 (N^3LO)  are provided by data
    //    member mBetaN.  The order of the expansion  mPertubativeOrder
    //    (defined as the 'n' in N^nLO) and the number of steps
    //    mIntegrationSteps  for the integration beyond LO are also kept
    //    as data member.
    //
    const double SXTH = 0.166666666666666;

    //
    //   Initial value, evolution distance and step size
    //
    double result = AS0;
    double LRRAT = log (R2/R20);
    double DLR = LRRAT / mIntegrationSteps;

    //
    //   Solution of the evolution equation depending on mPertubativeOrder
    //   (fourth-order Runge-Kutta beyond the leading order)
    //
    double XK0, XK1, XK2, XK3;
    if (mPertubativeOrder == 0) {
        result = AS0 / (1 + mBeta0[NF] * AS0 * LRRAT);
    }
    else if (mPertubativeOrder == 1) {
        for (int K1 = 1; K1 <= mIntegrationSteps; K1++) {
            XK0 = DLR * funBeta1 (result, NF);
            XK1 = DLR * funBeta1 (result + 0.5 * XK0, NF);
            XK2 = DLR * funBeta1 (result + 0.5 * XK1, NF);
            XK3 = DLR * funBeta1 (result + XK2, NF);
            result = result + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3);
        }
    }
    else if (mPertubativeOrder == 2) {
        for (int K1 = 1; K1 <= mIntegrationSteps; K1++) {
            XK0 = DLR * funBeta2 (result, NF);
            XK1 = DLR * funBeta2 (result + 0.5 * XK0, NF);
            XK2 = DLR * funBeta2 (result + 0.5 * XK1, NF);
            XK3 = DLR * funBeta2 (result + XK2, NF);
            result = result + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3);
        }
    }
    else if (mPertubativeOrder == 3) {
        for (int K1 = 1; K1 <= mIntegrationSteps; K1++) {
            XK0 = DLR * funBeta3 (result, NF);
            XK1 = DLR * funBeta3 (result + 0.5 * XK0, NF);
            XK2 = DLR * funBeta3 (result + 0.5 * XK1, NF);
            XK3 = DLR * funBeta3 (result + XK2, NF);
            result = result + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3);
        }
    }
    return result;
}

void AlphaStrong::betaFunction()
{
    //
    //   The function betaFunction() for the coefficients  mBeta0..mBeta3 of the
    //    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by
    //
    //        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ...
    //
    //    with  a_s = alpha_s/(4*pi).
    //
    //   The MSbar coefficients are written to the common-block BETACOM for
    //   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
    //
    //   The factors mCF, mCA and mTF  are data members.
    //    Beyond NLO the QCD colour factors are hard-wired in this routine,
    //    and the numerical coefficients are truncated to six digits.
    //

    const int NFMIN = 3;
    const int NFMAX = 6;

    //
    //  The full LO and NLO coefficients
    //
    double B00 =  11./3. * mCA;
    double B01 =  -4./3. * mTR;
    double B10 =  34./3. * mCA*mCA;
    double B11 = -20./3. * mCA*mTR - 4.* mCF*mTR;
    
    //
    //  Flavour-number loop and output to the array
    //
    for (int NF = NFMIN; NF <= NFMAX; NF++) {
        mBeta0[NF] = B00 + B01 * NF;
        mBeta1[NF] = B10 + B11 * NF;
        mBeta2[NF] = 1428.50 - 279.611 * NF + 6.01852 * NF*NF;
        mBeta3[NF] = 29243.0 - 6946.30 * NF + 405.089 * NF*NF + 1.49931 * NF*NF*NF;
    }
}

double AlphaStrong::rootFinder(double A0, double B0, double EPS, int MAXF, int MODE)
{
    //
    //   DZEROX taken from CERNLIB to find the zero of function findR0().
    //
    //     Based on
    //        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
    //        Guaranteed Convergence for Finding a Zero of a Function,
    //        ACM Trans. Math. Software 1 (1975) 330-345.
    //
    //        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
    //
    
    double result = 0;
    
    const double z1 = 1;
    const double half = z1/2;
    
    int IM1[3] = {0, 2, 3 };
    int IM2[3] = {0, -1, 3 };
    bool LMT[3];
    double C, D, FC, FD, FDA, FDB, P, Q, W;
    int IE;

    if (MODE != 1 && MODE != 2) {
        cout << "AlphaStrong::DZEROX(): Error, illegal mode (" << MODE << ")." << endl;
        return 0;
    }

    double FA = findR0(B0);
    double FB = findR0(A0);
    
    if (FA*FB > 0) {
        cout << "AlphaStrong::DZEROX(): Error, f(a) and f(b) have the same sign (a=" << A0 << ", b=" << B0 << ")." << endl;
        return 0;
   }

    double ATL = fabs(EPS);
    double B = A0;
    double A = B0;
    LMT[2] = true;
    int MF = 2;
L1:
    C = A;
    FC = FA;
L2:
    IE = 0;
L3:
    if (fabs(FC) < fabs(FB)) {
        if(C != A) {
            D=A;
            FD=FA;
        }
        A=B;
        B=C;
        C=A;
        FA=FB;
        FB=FC;
        FC=FA;
    }

    double TOL = ATL*(1+fabs(C));
    double H = half*(C+B);
    double HB = H-B;

    if (fabs(HB) > TOL) {
        if (IE > IM1[MODE]) {
            W = HB;
        }
        else {
            TOL = TOL*sign(z1,HB);
            P = (B-A)*FB;
            LMT[1] = (IE <= 1);
            if (LMT[MODE]) {
                Q = FA-FB;
                LMT[2] = false;
            }
            else {
                FDB = (FD-FB)/(D-B);
                FDA = (FD-FA)/(D-A);
                P = FDA*P;
                Q = FDB*FA-FDA*FB;
            }
            if (P < 0) {
                P = -P;
                Q = -Q;
            }
            if (IE == IM2[MODE]) P=P+P;
            if (P == 0 || P <= Q*TOL) {
                W = TOL;
            }
            else if (P < HB*Q) {
                W = P/Q;
            }
            else {
                W = HB;
            }
        }
        D = A;
        A = B;
        FD = FA;
        FA = FB;
        B = B+W;
        MF = MF+1;
        if (MF > MAXF) {
            cout << "AlphaStrong::DZEROX(): Too many function calls." << endl;
            return result;
        }
        FB = findR0(B);
        if (FB == 0 || sign(z1,FC) == sign(z1,FB)) goto L1;
        if (W == HB) goto L2;
        IE = IE+1;
        goto L3;
    }
    result = C;
    
    return result;
}
