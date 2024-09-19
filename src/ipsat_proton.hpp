/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015-2018
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


// If Fortran compiler is not available, uncommet this
// In that case, IPSAT12 does not work!
//#define USE_FORTRAN_IPSAT12

// In orde to use LCPT dipole (Dumitru, Mäntysaari, Paatelainen, 2021),
// One has to link uncomment this, and modify CMakeLists.txt files accordingly
//#define USE_LCPT_DIPOLE

#ifdef USE_LCPT_DIPOLE
#include <dipole_interpolation/dipoleamplitude.hpp>
#endif

// How are the hotspots distributed
enum Proton_shape
{
    GAUSSIAN,
    EXPONENTIAL,
};

// Quark structure
enum Structure
{
    QUARKS,     // Gaussians around quarks
    CENTER_TUBES,   // Quarks connected by flux tubes that merge at the center of the triangle
};

enum Fluctuation_shape
{
    FLUCTUATE_QUARKS,
};

struct IPsat_fit_parameteters
{
     // Default params: double C=2.2894; double mu0 = std::sqrt(1.1); double lambdag=0.08289; double Ag=2.1953; double mc=1.3528;
    double C;
    double mu0;
    double lambdag;
    double Ag;
    double mc;
    bool saturation;
};

class Ipsat_Proton : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    double Amplitude(double xpom, Vec q1, Vec q2);  // this is also in DipoleAmplitude, a bit overlap...
    
    double Amplitude_Tp(double xpom, double r, double tp);  // Dipole ampliutde at fixed tp
    
    double Amplitude_bint(double xpom, double r);
    double Amplitude_sqr_bint(double xpom, double r);
    
    double xg(double x, double r);
    double LogDerivative_xg(double x, double r); // d ln xg / d ln (1/x)
    
    // Setup the target. In practice sample nucleon positions from Woods Saxon
    void InitializeTarget();
    
    Ipsat_Proton(DGLAPDist *gd);    // Use global gdist
    Ipsat_Proton(Ipsat_version version);
    Ipsat_Proton(Ipsat_version version, IPsat_fit_parameteters params);
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
    
    void SetQsFluctuation(double s);    // Set sigma for ln Q_s fluctuation
    double GetQuarkQsFluctuation(unsigned int i);   // Get Q_s fluctuation for the given quark
    
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
    
#ifdef USE_LCPT_DIPOLE
    LCPT_Dipole *lcpt_dipole;
#endif
    
    bool skewedness;    // Enable skewedness in dipole amplitude, multiplies xg
    
    bool origin_at_center_of_mass;    // If true, move quarks s.t. their center of mass is at b=0
    
    gsl_integration_workspace *intworkspace_zint;   // Used to calculate z int of exponential/fluxtube distribution
};

#endif /* ipsat_proton_hpp */
