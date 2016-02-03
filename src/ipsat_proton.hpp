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

enum Proton_shape
{
    GAUSSIAN,
    EXPONENTIAL
};

enum Ipsat_version
{
    IPSAT06,    // KMW hep-ph/0606272
    IPSAT12     // Rezaeian et al, 1212.2974
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
    ~Ipsat_Proton();
    
    std::vector<Vec>& GetQuarks();
    std::vector<double> GetRadii();
    
    // Set size parameters, affect to next InitializeTarget() call
    void SetProtonWidth(double bp_);
    void SetQuarkWidth(double bq_);
    
    void SetShape(Proton_shape s);
    
    std::string InfoStr();
    
    void SetSkewedness(bool s);     // Set whether we include skewedness when calculating amplitude
    
    void SampleQsFluctuations();
    
    double GetQsFluctuation(double x, double y);   // return Exp(f(x,y)) that multiplies xg
    void SetQsFluctuation(double s);    // Set sigma for ln Q_s fluctuation
    
    void SetFluctuationShape(Fluctuation_shape s);
    Fluctuation_shape GetFluctuationShape();
    
private:
    void Init();
    
    double Skewedness(double lambda);
    DGLAPDist *gdist;    // DGLAP evolved xg
    // gdist.Gluedist() returns Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
    std::vector<Vec> quarks;    // Quark positions
    std::vector<double> quark_bp;
    double maxr;
    double B_p;     // Central value for the proton gaussian width
    double B_q;     // Central value for the quark gaussian width
    bool saturation;
    
    
    
    double QuarkThickness(double r, int i); // Quark density profile for quark i, distance r from its origin
    
    double GaussianRadiusDistribution(double r);    // Used to sample Gaussian radius
    double ExponentialDistribution(double x, double y, double z);   // Exponential distribution for a vector
    
    bool allocated_gdist;   // True if we have allocated memory for gdist in here
    
    Proton_shape shape;
    
    // Q_s luctuations
    Fluctuation_shape fluctuation_shape;
    double Qs_fluctuation_sigma;        // Width of the ln Q_s^2/<Q_s^2> distribution
    std::vector< std::vector< double > > qs_fluctuation;    // Grid of Q_s fluctuations
    std::vector< double > qs_fluctuation_coordinates;       // Grid points for Q_s^2 fluctuation
    std::vector< double > qs_fluctuations_quarks;           // Q_s^2 fluctuations for each quark
    
    
    Ipsat_version ipsat;
    
    bool skewedness;    // Enable skewedness in dipole amplitude, multiplies xg
};

#endif /* ipsat_proton_hpp */
