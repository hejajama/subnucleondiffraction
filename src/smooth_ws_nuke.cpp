/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for dipole-smooth nucleus scattering for testing/comparisons
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "smooth_ws_nuke.hpp"
#include "mz_ipsat/dipoleamplitude.hpp"
#include <tools/tools.hpp>
#include <tools/interpolation.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

// IPsat 2012
extern "C" {
    double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};

int IPSAT12_NUKE_PAR = 2;    // m_c=1.4 GeV

using Amplitude::SQR;

const double NC=3.0;

using namespace std;

Smooth_ws_nuke::Smooth_ws_nuke(int A_)
{
    mzipsat = new MZ_ipsat::DipoleAmplitude(2.146034445992, 1.1, 0.09665075464199, 2.103826220003, 1.351650642298);
    mzipsat->SetSaturation(true);
    
    A=A_;
    InitializeWSDistribution(A);
    
    // Initialize interpolator
    vector<double> bvals;
    vector<double> tavals;
    for (double b=0; b<100; b+=0.1)
    {
        bvals.push_back(b);
        tavals.push_back(T_A(b, A));
    }
    T_A_interpolator = new Interpolator(bvals, tavals);
    T_A_interpolator->SetOverflow(0);
    T_A_interpolator->SetUnderflow(0);
    T_A_interpolator->SetFreeze(true);
}

Smooth_ws_nuke::~Smooth_ws_nuke()
{
    delete T_A_interpolator;
}
double Smooth_ws_nuke::Amplitude(double xpom, double q1[2], double q2[2] )
{
        
    // Very crude toy model: neglet different positions for quarks, only take impact parameter dependence from the Woods Saxon
    // Basically KT hep-ph/0304189 Eq. (41)
    
    double r = std::sqrt( SQR(q1[0]-q2[0]) + SQR(q1[1]-q2[1]) );
    
    // Take the nucleon density at the geometric center of the two quarks
    double b = std::sqrt( SQR( (q1[0]+q2[0])/2.0 ) + SQR( (q1[1]+q2[1])/2.0) );
    
    return 1.0 - std::exp( -r*r * M_PI*M_PI / (2.0 * NC) * mzipsat->Alphas_xg(xpom, mzipsat->MuSqr(r)) * A * T_A_interpolator->Evaluate(b));
    
    
   /* if (ipsat == IPSAT06)
    {
        return 1.0 - exp( -A*gdist.Gluedist(xpom, r*r) *r*r* T_A_interpolator->Evaluate(b));
    }
    else if (ipsat == IPSAT12)
    {*/
    /*
        // dipole_amplitude(xBj, r, b, parameterSet) gives amplitude 2(1 - exp(c*T(p)))
        // We have to calculate "gluedist" as in case of ipsat06
        double tmpb=0;
        // par 1: m_c=1.27,   2: m_c=1.4
        double n = dipole_amplitude_(&xpom, &r, &tmpb, &IPSAT12_NUKE_PAR)/2.0;
        
        
        double c = std::log(1.0-n);
        
        if (std::isnan(c) or std::isinf(c))
        {
            // We have so large r, that basically n=1 and c blows up, these should not matter
            // as wave function cuts these out anyway, but we can just set amplitude to 1
            return 1.0;
        }
        
        double tp = 1.0/(2.0*M_PI*4.0)*std::exp(- tmpb*tmpb / (2.0*4.0));
        c /= tp;
        
        /*double skew=1.0; // now c contains xg that is modified by skewedness correction if enabled
        if (skewedness and tmpr < MAXR_SKEW)
        {
            double skew_lambda = LogDerivative_xg(xpom, r.Len());
            if (std::isnan(skew_lambda))
            {
                std::cerr << "skew nan at xpom=" << xpom << " r = " << r.Len() << std::endl;
                skew=1.0;
            }
            else
            skew = Skewedness(skew_lambda);
        }
    
        return 1.0 - std::exp( A *T_A_interpolator->Evaluate(b) * c ); // c contains the - sign
    //}

    */
    
    
    // As T_A is normalized to unity, we get no extra normalizatino factor for it
    
}

std::string Smooth_ws_nuke::InfoStr()
{
    std::stringstream ss;
    ss << "#Optigal Glauber nucleus, A=" << A << endl;;
    return ss.str();
}
