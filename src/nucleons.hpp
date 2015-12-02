/*
 * Nucleus that consists of nucleons
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef nucleons_hpp
#define nucleons_hpp

#include "dipole.hpp"
#include <vector>
#include "ipsat_proton.hpp"
#include "vector.hpp"
#include "gdist_dglap.hpp"

class Nucleons : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    void InitializeTarget();
    
    // Set size parameters, affect to next InitializeTarget() call
    void SetProtonWidth(double bp_);
    void SetQuarkWidth(double bq_);
    
    Nucleons();
    double WS_unnorm(double r );
    
private:
    int A;
    double B_p;       // Proton size parameter
    double B_q;     // Quark size parameter
    std::vector<Ipsat_Proton> nucleons; // Todo: change to support general proton/nucleon class?
    std::vector<Vec> nucleon_positions;
    
    double ws_delta;
    double ws_ra;
    
    DGLAPDist gdist;    // Use when create nucleus from ipsat protons, so that not every proton has to
    // load the gluon distribution
    
};

#endif /* nucleons_hpp */
