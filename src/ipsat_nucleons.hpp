/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a nucleus that consist of nucleons described ty IPsat
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef ipsat_nucleons_hpp
#define ipsat_nucleons_hpp

#include "dipole.hpp"
#include <vector>
#include "vector.hpp"
#include "gdist_dglap.hpp"

class Ipsat_Nucleons : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    // Setup the target. In practice sample nucleon positions from Woods Saxon
    void InitializeTarget();
    
    std::vector<Vec>& GetNucleons();
    
    Ipsat_Nucleons();
    
    double MaxR();  // Max radius used to sample nucleons, depends on A, in GeV^-1
    
    double Tp(double b);   // Proton thickness function
    
    double WS_unnorm(double r);  // WS distribution normalized to unity
    
    void SetSaturation(double s);
    
private:
    int A;
    double mindist; // How close nucleons can be
    double ws_delta;   // Delta parameter in the WS distribution, units GeV^-1
    double ws_ra;       // Nuclear radius in WS distribution
    
    double B_p;     // Parameters proton transverse size, goes into the b dependence if ipsat
    
    DGLAPDist gdist;    // DGLAP evolved xg
    // gdist.Gluedist() returns Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
    
    bool saturation;    // Turn saturation on/off. Without saturation use linearized dipole amplitude ~r^2
    
    std::vector<Vec> nucleons; // Positions of the centers of the nucleons
};


#endif /* ipsat_nucleons_hpp */
