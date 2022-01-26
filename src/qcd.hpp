#ifndef _QCD_HPP
#define _QCD_HPP

 
enum RunningAlphas
{
    RUNNING,
    FIXED
};

/**
  * Virtual photon polarizatoin states
  */
enum Polarization
{
    L,  ///<< Longitudinal polarization
    T   ///<< Transverse polarization
};

enum Parton
{
    UVAL,   ///< valence u
    DVAL,   ///< valence d
    USEA,   ///< sea u
    DSEA,   ///< sea d
    U,      ///< u (valence + sea)
    D,      ///< d (valence + sea)
    S,      ///< strange
    C,      ///< charm
    B,      ///< bottom
    G,      ///< gluon
    UBAR,   ///< anti-u
    DBAR,   ///< anti-d
    SBAR,   ///< anti-s
    CBAR,   ///< anti-c
    BBAR,   ///< anti-b
    LIGHT   ///< light quarks = u+d+s
};

const double ALPHA_e = 1.0/137.035999679;

/**
 * QCD functions and constants
 *
 * Strong coupling constant, Nc, Nf etc. parameters
 */

class QCD
{
    public:
        QCD();
        
        /**
         * Strong coupling constant
         * @param musqr scale [GeV^2]
         */
        double Alphas(double musqr);
        /**
         * Number of quark flavors
         */
        int Nf();
        /**
         * Number of colors
         */
        int Nc();

        /**
         * Fundametal Casimir
         */
        double Cf();

        /**
         * Set running/fixed coupling
         */
        void SetRunningCoupling(RunningAlphas rc_);
    private:
        int nc;
        int nf;
        double lqcd;    // lambda_QCD
        double maxalpha;    // freeze coupling to maxalpha in infrared
        RunningAlphas rc;
        

};

double T_A(double b, int A);    // Transverse density profile of the nucleus normalized by A
                                                                // Unit: [b]=1/GeV!
void InitializeWSDistribution(int A);


#endif

