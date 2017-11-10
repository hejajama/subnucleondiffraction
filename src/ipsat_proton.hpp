/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */


#ifndef ipsat_proton_hpp
#define ipsat_proton_hpp

#include <gsl/gsl_integration.h>

#include "dipole.hpp"
#include <vector>
#include <string>
#include "vector.hpp"
#include "gdist_dglap.hpp"
#include "vector.hpp"
#include "mz_ipsat/dipoleamplitude.hpp"

// How are the hotspots distributed
enum Proton_shape
{
    GAUSSIAN,
    EXPONENTIAL,
    ALBACETE,   // implement 1605.09176
};

// Quark structure
enum Structure
{
    QUARKS,     // Gaussians around quarks
    CENTER_TUBES,   // Quarks connected by flux tubes that merge at the center of the triangle
};

enum Ipsat_version
{
    IPSAT06,    // KMW hep-ph/0606272
    IPSAT12,     // Rezaeian et al, 1212.2974
    MZ
};

enum Fluctuation_shape
{
    FLUCTUATE_QUARKS,
    LOCAL_FLUCTUATIONS
};

class Ipsat_Proton : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    double Amplitude(double xpom, Vec q1, Vec q2);  // this is also in DipoleAmplitude, a bit overlap...
    
    double xg(double x, double r);
    double LogDerivative_xg(double x, double r); // d ln xg / d ln (1/x)
    
    // Setup the target. In practice sample nucleon positions from Woods Saxon
    void InitializeTarget();
    
    Ipsat_Proton();
    Ipsat_Proton(DGLAPDist *gd);    // Use global gdist
    Ipsat_Proton(Ipsat_version version);
    ~Ipsat_Proton();
    
    std::vector<Vec>& GetQuarks();
    std::vector<double> GetRadii();
    
    // Set size parameters, affect to next InitializeTarget() call
    void SetProtonWidth(double bp_);
    void SetQuarkWidth(double bq_);
    
    void SetShape(Proton_shape s);
    void SetStructure(Structure s);
    
    std::string InfoStr();
    
    void SetSkewedness(bool s);     // Set whether we include skewedness when calculating amplitude
    
    void SampleQsFluctuations();
    
    double GetQsFluctuation(double x, double y);   // return Exp(f(x,y)) that multiplies xg
    void SetQsFluctuation(double s);    // Set sigma for ln Q_s fluctuation
    
    void SetFluctuationShape(Fluctuation_shape s);
    Fluctuation_shape GetFluctuationShape();
    
    double FluxTubeThickness(Vec b); // Density in the flux tube model
    
    double QuarkThickness(double r, int i); // Quark density profile for quark i, distance r from its origin
    
    void SetFluxTubeNormalization(double n);
    
    double Density(Vec b);    // Density (T_b) at given point b
    
    void SetQuarkCenterOfMassToOrigin(bool s);
    
private:
    void Init();
    
    int number_of_quarks;
    
    double Skewedness(double lambda);
    DGLAPDist *gdist;    // DGLAP evolved xg
    // gdist.Gluedist() returns Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
    std::vector<Vec> quarks;    // Quark positions
    std::vector<Vec> quarks3d;  // Quark positions including z coordinate
    Vec center;                 // Fermat point of the quark triangle projected to z=0
    Vec center3d;               // Fermat poitn of the quark triangle in 3d
    std::vector<double> quark_bp;
    double maxr;
    double B_p;     // Central value for the proton gaussian width
    double B_q;     // Central value for the quark gaussian width
    bool saturation;
    
    
    double fluxtube_normalization; // Normalization factor for FluxTube
    void NormalizeFluxTubeThickness();  // Calculate normalization factor, do this after quark positions
    // are sampled
    
    
    
    double GaussianRadiusDistribution(double r);    // Used to sample Gaussian radius
    double ExponentialDistribution(double x, double y, double z);   // Exponential distribution for a vector
    double RepulsiveGaussianDistribution(std::vector<Vec> quarks, double rc); // Distribution from 1605.09176
    
    bool allocated_gdist;   // True if we have allocated memory for gdist in here
    
    Proton_shape shape;
    
    Structure proton_structure;
    
    // Q_s luctuations
    Fluctuation_shape fluctuation_shape;
    double Qs_fluctuation_sigma;        // Width of the ln Q_s^2/<Q_s^2> distribution
    std::vector< std::vector< double > > qs_fluctuation;    // Grid of Q_s fluctuations
    std::vector< double > qs_fluctuation_coordinates;       // Grid points for Q_s^2 fluctuation
    std::vector< double > qs_fluctuations_quarks;           // Q_s^2 fluctuations for each quark
    
    
    Ipsat_version ipsat;
    
    MZ_ipsat::DipoleAmplitude *mzipsat;
    
    bool skewedness;    // Enable skewedness in dipole amplitude, multiplies xg
    
    bool origin_at_center_of_mass;    // If true, move quarks s.t. their center of mass is at b=0
    
    gsl_integration_workspace *intworkspace_zint;   // Used to calculate z int of exponential/fluxtube distribution
};

#endif /* ipsat_proton_hpp */
