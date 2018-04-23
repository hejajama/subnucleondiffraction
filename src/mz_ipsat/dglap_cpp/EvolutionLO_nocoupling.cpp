//==============================================================================
//  EvolutionLO.cpp
//
//==============================================================================
#include "EvolutionLO_nocoupling.h"
#include <cmath>
#include <cstdlib>

#define PR(x) cout << #x << " = " << (x) << endl;

inline complex<double> SQR(complex<double> x)
{
    return x*x;
}

double EvolutionLO_gluon::alphasxG(double x, double Q2, double mu0, int coupling, double Ag,
                             double lambdag, double As, double lambdas)
{
    double val = xG(x, Q2, mu0, coupling, Ag, lambdag, As, lambdas);
    return val*mAlphaStrong->value(sqrt(Q2));
    
}

double EvolutionLO_gluon::xG(double x, double Q2, double mu0, int coupling, double Ag,
                       double lambdag, double As, double lambdas)
{
    if (mUseLookupTable && !mLookupTableIsFilled) {
        cerr << "EvolutionLO_gluon::xG(): Warning, use of lookup table requested\n"
        << "                           but table is not setup. First use \n"
        << "                           generateLookupTable() to fill table.\n"
        << "                           Will fall back to numeric calculation."
        << endl;
    }
    if (mUseLookupTable and x > mTableMinX and x < mTableMaxX and Q2 > mTableMinQ2 and Q2 < mTableMaxQ2)
        return xG_Interpolator(x, Q2);
    
    const double FourPi = 4.*M_PI;
    
    mMUR = mu0;               // input mu_r in GeV
    
    //
    //   These are provided by AlphaStrong ensuring
    //   consistency between evolution and alpha_s.
    //
    if (mAlphaStrong->order() != 0) {
        cout << "EvolutionLO::xG(): Fatal error, alpha_s is in order "
             << mAlphaStrong->order()
             << " but EvolutionLO only support order = 0. Stop here." << endl;
        exit(1);
    }
    mFR2 = mAlphaStrong->ratioFR2();
    mASMUR = mAlphaStrong->value(mu0);
    mMC = mAlphaStrong->massCharm();
    mMB = mAlphaStrong->massBottom();
    mMT = mAlphaStrong->massTop();

    mALPS = mAlphaStrong->value(mMUR)/FourPi;
    mALPC = mAlphaStrong->value(mMC)/FourPi;
    mALPB = mAlphaStrong->value(mMB)/FourPi;
    mALPT = mAlphaStrong->value(mMT)/FourPi;
    double alpq = mAlphaStrong->value(sqrt(Q2))/FourPi;
    
    //
    //   Q2 and x dependent quantities
    //
    double ax = log(x);
    double ex = exp(-mC * ax);
    
    //
    //  integration length parameter for the mellin inversion
    //
    const int nmax = 136;

    //
    //   Gluon density and output
    //
    complex<double> fn[137];
    reno(fn, alpq, nmax, coupling, Ag, lambdag, As, lambdas);
    
    long double fun = 0;   // long double is needed
    long double fz;
    complex<double> xnm,cex;
    for (int i=1; i <= nmax; i++) {
        xnm = (mC - mN[i]+1.) * ax;
        cex = exp(xnm) / M_PI * mCC;
        fz = imag(fn[i]*cex);
        fun = fun + mWN[i] * fz;
    }
    double pa = fun * ex;
    
    return pa;
}

void EvolutionLO_gluon::reno(complex<double> *fn, double alpq, int nmax, int coupling,
                       double ag, double lambdag, double as, double lambdas)
{
    if (coupling==1)
    {
        cerr << "Using EvolutionLO_gluon with coupling=1, this can not be right!" << endl;
        exit(1);
    }
    //
    //   Mellin-n space Q**2 - evolution of the gluon at LO
    //
    //    The moments are calculated on an array of moments,   mN,  suitable
    //    for a (non-adaptive) Gauss quadrature.
    //
    //    Currently this takes the simplest possible fit form:
    //    xg = A_g x^(-lambdag) (1-x)^(5.6), following Amir&Raju
    //
    for (int k1 = 1; k1 <= nmax; k1++) {
        //
        //   Input moments of the parton densities
        //   at the low scale
        //
        complex<double> xn =  mN[k1];
        complex<double> zero;
        complex<double> one(1.,0.);
        complex<double> two(2.,0.);
        
        // x*g = A_g x^(-lambda_g) (1-x)^5.6
        // -one and +one as we parametrize g, not x*g
        //complex<double> gln = ag * beta(xn - lambdag - one, 5.6+one);
        
        // x*g = A_g x^(-lambda_g) (1-x)^6
        // Mellin transformed g analytically
        
        
        
        complex<double> gln = ag * (
            1.0 / (xn + 5.0 - lambdag)
            - 6.0 / (xn + 4.0 - lambdag)
            + 15.0 / (xn + 3.0 - lambdag)
            - 20.0 / (xn + 2.0 -lambdag)
            + 15.0 / (xn + 1.0 - lambdag)
            - 6.0 / (xn - lambdag)
            + 1.0 / (xn - lambdag - 1.0)
                                    );
       
	// Test: x*g = A_g x^(-lambda_g) (1-x)^1
	 //complex<double> gln = ag * ( 1.0 / ( (-1.0 + xn - lambdag) * (xn - lambdag) ) ); 

	//(1-x)^2
	//complex<double> gln = ag * 2.0 / ( (-1.0 + xn - lambdag)*(xn - lambdag) * (1.0 + xn - lambdag) );
  
        
        int f;
        double xl, s, alp;
        complex<double> ep, gl;
        
        if (alpq >= mALPC) {   // evolution below the charm threshold
            f = 3;
            xl = mALPS / alpq;
            
            s = log(xl);
            
            ep = exp(-mAP[k1][f]*s);
            
            gl = gln;
            
           
            gln = ep * gl;
           
            
        }
        else if ((alpq < mALPC) && (alpq >= mALPB)) {  // between thresholds
            f = 3;
            xl = mALPS / mALPC;
            
            s   = log(xl);

            ep  = exp(-mAP[k1][f]*s);
            
            gl = gln;
            
            
            gln = ep * gl;
            
            
            
            f = 4;
            xl = mALPC / alpq;
            
            
            s   =  log(xl);
            ep  = exp(-mAP[k1][f]*s);
            
            gl = gln;
            
            
            gln = ep * gl;
            
        }
        else if (alpq < mALPB) {    // above bottom threshold
            f = 3;
            xl = mALPS / mALPC;
            
            s   = log (xl);
            ep  = exp (-  mAP[k1][f]*s);
    
            gl = gln;
            
            
            gln = ep * gl;
            
            f = 4;
            alp = mALPB;
            xl = mALPC / mALPB;
            
            s   = log (xl);
           
            ep  = exp(-mAP[k1][f]*s);
            
            gl = gln;
            
            
            gln = ep * gl;
            
            f = 5;
            xl = mALPB / alpq;
            
            
            s  = log(xl);
            
            ep  = exp(-mAP[k1][f]*s);
            
            gl = gln;
            
            
            gln = ep * gl;
        }
        fn[k1] = gln;
    }
}


EvolutionLO_gluon::EvolutionLO_gluon(AlphaStrong* alphas)
{
    //
    //  Get alpha_s
    //
    if (!alphas) {
        cout << "EvolutionLO::EvolutionLO(): Fatal error. Need valid AlphaStrong object. Stop here." << endl;
        exit(1);
    }
    mAlphaStrong = alphas;

    
    // Interpolation options
    mTableMinQ2 = 0.6;
    mTableMaxQ2 = 1e8;
    mTableMinX = 1e-8;
    mTableMaxX = 0.05;
    mLookupTableIsFilled = false;
    mUseLookupTable = false;
    
    //
    //   Initialization of support points in n-space and weights for the
    //   Gauss quadrature and of the anomalous dimensions for the RG
    //   evolution at these n-values.
    //
    
    //
    //  Weights and support points for nomalized 8 point gauss quadrature
    //
    double wz[9] = {0, 0.101228536290376,0.222381034453374,0.313706645877887,
        0.362683783378362,0.362683783378362,0.313706645877887,
        0.222381034453374,0.101228536290376};
    double zs[9] = {0, -0.960289856497536,-0.796666477413627,-0.525532409916329,
        -0.183434642495650,0.183434642495650,0.525532409916329,
        0.796666477413627,0.960289856497536};
    //
    //  Integration contour parameters
    //
    double down[18] = {0, 0., 0.5, 1., 2., 3., 4., 6., 8.,
        1.e1, 1.2e1, 1.5e1, 1.8e1, 2.1e1, 2.4e1, 2.7e1, 3.e1, 3.3e1};
    double up[18];
    mC = 0.8;
    double phi = M_PI * 3./4.;
    double co = cos(phi);
    double si = sin(phi);
    mCC = complex<double>(co, si);
    for (int i=1; i <=16; i++) up[i] = down[i+1];
    up[17] = 36.;
    
    //
    //  Support points and weights for the gauss integration
    //  (the factor (up-down)/2 is included in the weights)
    //
    int k = 0;
    double sum, diff, z;
    for (int i=1; i <=17; i++) {
        sum  = up[i] + down[i];
        diff = up[i] - down[i];
        for (int j=1; j <=8; j++) {
            k++;
            z = 0.5 * (sum + diff * zs[j]);
            mWN[k] = diff / 2.* wz[j];
            mN[k]  = complex<double>(mC+co*z+1.,si*z);
        }
    }
    anom();
 }

EvolutionLO_gluon::~EvolutionLO_gluon()
{
    if (mAlphaStrong) delete mAlphaStrong;
    if (mLookupTableIsFilled) {
        for (unsigned int i = 0; i < mNumberOfNodesInQ2; ++i)
        delete [] mLookupTable[i];
        delete [] mLookupTable;
    }
}




void EvolutionLO_gluon::generateLookupTable(double mu0, int coupling, double Ag,
                                            double lambdag, double As, double lambdas, int nx, int nq2 )
{
    //
    //  Delete old lookup table (if it exists)
    //
    if (mLookupTableIsFilled) {
        for (unsigned int i = 0; i < mNumberOfNodesInQ2; ++i)
         delete [] mLookupTable[i];
        delete [] mLookupTable;
    }
    
    mNumberOfNodesInX = nx;
    mNumberOfNodesInQ2 = nq2;
    
    //
    //  Create new table
    //
    mLookupTable = new double*[mNumberOfNodesInQ2];
    for(unsigned int i = 0; i < mNumberOfNodesInQ2; ++i)
    mLookupTable[i] = new double[mNumberOfNodesInX];
    
    //
    //  Fill lookup table
    //
    int printEach = mNumberOfNodesInQ2*mNumberOfNodesInX/10;
    int nCount = 0;
//    cout << "#DglapEvolution: generating lookup table "; cout.flush();
    for (unsigned int i = 0; i < mNumberOfNodesInQ2; i++) {
        double Q2 = pow(mTableMaxQ2/mTableMinQ2,i*1./(mNumberOfNodesInQ2-1.))*mTableMinQ2;
        for (unsigned int j = 0; j < mNumberOfNodesInX; j++) {
            double x = pow(mTableMaxX/mTableMinX,j*1./(mNumberOfNodesInX-1.))*mTableMinX;
            mLookupTable[i][j] = xG(x, Q2, mu0,  coupling,  Ag,
                                     lambdag,  As,  lambdas);
            nCount++;
  //          if (nCount%printEach == 0) cout << '.'; cout.flush();
        }
    }
 //   cout << " done" << endl;
    mLookupTableIsFilled = true;
}

double EvolutionLO_gluon::xG_Interpolator(double x, double Q2)
{
    int q2steps = mNumberOfNodesInQ2-1;
    int xsteps = mNumberOfNodesInX-1;
    
    //
    //  Position in the Q2 grid
    //
    double realq = q2steps * log(Q2/mTableMinQ2)/log(mTableMaxQ2/mTableMinQ2);
    
    int n_q2 = static_cast<int>(realq);
    if (n_q2 <= 0) {n_q2 = 1;}
    if (n_q2 > q2steps-2) {n_q2 = n_q2-2;}
    
    //
    //  Position in the x grid
    //
    double realx = xsteps * log(x/mTableMinX)/log(mTableMaxX/mTableMinX);
    
    int n_x = static_cast<int>(realx);
    if (n_x <= 0) { n_x = 1;}
    if (n_x > xsteps-2){ n_x = n_x-2;}
    
    //
    //  Starting the interpolation
    //
    double arg[4];
    for (int i=0; i < 4; i++) { arg[i] = n_x-1+i;}
    
    double fu[4], fg[4];
    for (int i=0; i < 4; i++) {
        fu[0] = mLookupTable[n_q2-1+i][n_x-1];
        fu[1] = mLookupTable[n_q2-1+i][n_x];
        fu[2] = mLookupTable[n_q2-1+i][n_x+1];
        fu[3] = mLookupTable[n_q2-1+i][n_x+2];
        fg[i] = luovi(fu,arg,realx);
    }
    for (int i=0; i < 4; i++) { arg[i] = n_q2-1+i;}
    
    return luovi(fg, arg, realq);
}

double EvolutionLO_gluon::luovi(double f[4], double arg[4], double z)
{
    double cof[4];
    for (int i=0; i < 4; i++ ) cof[i]=f[i];
    
    for (int i=0; i < 3 ; i++) {
        for (int j=i; j < 3 ; j++) {
            int jndex = 2 - j;
            int index = jndex + 1 + i;
            cof[index]= (cof[index]-cof[index-1])/(arg[index]-arg[jndex]);
        }
    }
    
    double sum=cof[3];
    
    for (int i=0; i < 3; i++ ) {
        int index = 2 - i;
        sum = (z-arg[index])*sum + cof[index];
    }
    
    return sum;
}



void EvolutionLO_gluon::anom() {
    //
    //   Anomalous dimensions for leading order evolution of parton densities.
    //   The moments are calculated on an externally given array of mellin
    //    moments,   mN, suitable for a (non-adaptive) quadrature.
    //
    //   Present version: the number of active flavours in the factorization
    //    is fixed to ff=3, in the beta-function it varies between f=3 and
    //    f=5. The dimension of the moment array is 136.
    //

    double b0,  b02;
    complex<double>  ggi, ggf;
    complex<double> xn, xap;
    complex<double>  gg;
    
    for (int k1=1; k1 <= 136; k1++) {
        xn = mN[k1];
        anCalc(ggi, ggf, xn);
        for (int k2=3; k2 <= 5; k2++) {
            double f = k2;
            //  anomalous dimensions and related quantities in leading order
            b0 = 11.- 2./3.* f;
            b02 = 2.* b0;
            gg = ggi + f * ggf;
            xap = gg / b02;
            mAP[k1][k2] = xap;
        }
    }
}

void EvolutionLO_gluon::anCalc( complex<double>& ggi,
              complex<double>& ggf, complex<double>& xn)
{
    complex<double> xns = xn * xn;
    complex<double> xn1 = xn + 1.;
    complex<double> xn2 = xn + 2.;
    complex<double> xnm = xn - 1.;

    //
    //  Leading order
    //
    complex<double> cpsi = psiFunction(xn1) + 0.577216;
    ggi = -22.- 24./(xn * xnm) - 24./(xn1 * xn2) + 24.* cpsi;
    ggf = 4./3.;
}

complex<double> EvolutionLO_gluon::psiFunction(complex<double> z)
{
    //
    //  psi - function for complex argument
    //
    complex<double> sub;

    while (real(z) < 10) {
        sub = sub - 1./z;
        z = z + 1.;
    }
    complex<double> rz = 1./z;
    complex<double> dz = rz * rz;
    complex<double> result = sub + log(z) - rz/2.- dz/2520. * ( 210.+ dz * (-21.+10.*dz ));
    return result;
}

complex<double> EvolutionLO_gluon::lngam(complex<double> x)
{
    complex<double> result = (x - complex<double>(0.5,0.0)) * log(x) - x +
    complex<double>(0.91893853,0.0) + complex<double>(1.0,0.0)/(complex<double>(12.0,0.0)* x)
    * (complex<double>(1.0,0.0)- complex<double>(1.0,0.0)/(complex<double>(30.,0.)* x*x)
       * (complex<double>(1.,0.)- complex<double>(1.,0.0)/(complex<double>(3.5,0.0) * x*x)
          * (complex<double>(1.,0.0)- complex<double>(4.,0.0)/(complex<double>(3.,0.0)* x*x))));
    
    return result;
}

complex<double> EvolutionLO_gluon::beta(complex<double> z1, complex<double> z2)
{
    complex<double> sub;
    
    while (real(z1) < 15) {
        sub = sub + log((z1+z2) / z1);
        z1 = z1 + 1.;
    }
    
    while (real(z2) < 15) {
        sub = sub + log ((z1+z2) / z2);
        z2 = z2 + 1.;
    }
    
    complex<double> lg1 = lngam(z1);
    complex<double> lg2 = lngam(z2);
    complex<double> lg12 = lngam(z1 + z2);
    return exp(lg1 + lg2 - lg12 + sub);
}
