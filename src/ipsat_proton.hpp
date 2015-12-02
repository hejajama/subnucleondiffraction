/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */


#ifndef ipsat_proton_hpp
#define ipsat_proton_hpp


#include "dipole.hpp"
#include <vector>
#include <string>
#include "vector.hpp"
#include "gdist_dglap.hpp"

class Ipsat_Proton : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    double Amplitude(double xpom, Vec q1, Vec q2);  // this is also in DipoleAmplitude, a bit overlap...
    
    
    // Setup the target. In practice sample nucleon positions from Woods Saxon
    void InitializeTarget();
    
    Ipsat_Proton();
    Ipsat_Proton(DGLAPDist *gd);    // Use global gdist
    ~Ipsat_Proton();
    
    std::vector<Vec>& GetQuarks();
    std::vector<double> GetRadii();
    
    // Set size parameters, affect to next InitializeTarget() call
    void SetProtonWidth(double bp_);
    void SetQuarkWidth(double bq_);
    
    std::string InfoStr();
    
    
private:
    DGLAPDist *gdist;    // DGLAP evolved xg
    // gdist.Gluedist() returns Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
    std::vector<Vec> quarks;    // Quark positions
    std::vector<double> quark_bp;
    double maxr;
    double B_p;     // Central value for the proton gaussian width
    double B_q;     // Central value for the quark gaussian width
    bool saturation;
    
    double QuarkThickness(double r, int i); // Quark density profile for quark i, distance r from its origin
    
    double RadiusDistribution(double r); // Distribution used to sample distances of quarks from the origin
    
    bool allocated_gdist;   // True if we have allocated memory for gdist in here
};

#endif /* ipsat_proton_hpp */
