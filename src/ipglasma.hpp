/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a IPglasma nucleus
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef ipglasma_hpp
#define ipglasma_hpp

#include <string>
#include <vector>
#include "wilsonline.hpp"
#include "dipole.hpp"

enum WilsonLineDataFileType
{
    BINARY,
    TEXT
};

class IPGlasma : public DipoleAmplitude {
public:
    
    IPGlasma(std::string fname);
    // Specify the step size, note that here step is in fm!
    IPGlasma(std::string fname, double step, WilsonLineDataFileType=TEXT);
    
    int LoadData(std::string fname, double step, WilsonLineDataFileType type=TEXT);
    int LoadBinaryData(std::string fname, double step);
    
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2]);
    double AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] );
    double Amplitude(double xpom, Vec q1, Vec q2);

    WilsonLine& GetWilsonLine( double x, double y); // Find Wilson line that corresponds to the coordinate
    
    std::string InfoStr();
    
    double MinX();
    double MaxX();
    
    double XStep(); // Grid spacing in x

    void SetPeriodicBoundaryConditions(bool s) { periodic_boundary_conditions = s; }
    
    std::vector<double> &GetXCoordinates();

    void SetSchwinger(bool s, double rc=0);
    void ApplyPeriodicBoundaryConditions(double q[2]); 

    void SetPeripheralExclusion(double B, double R_A, double th) { AA_impact_parameter=B; AA_nuclear_radius=R_A; AA_theta = th; }

private:
    
    
    std::vector< double > xcoords;
    std::vector< double > ycoords;
    std::vector< WilsonLine  > wilsonlines;

    bool schwinger; // true if we use schwinger mechanism
    double schwinger_rc; // Use Schwinger for dipoles larger than schwinger_rc

    bool periodic_boundary_conditions;
    double AA_impact_parameter; // If >0, consider peripheral collisions and remove points that are closer to the other nucleus than AA_nuclear_radius.
    double AA_nuclear_radius;    
    double AA_theta; // Azimuthal angle of the other nucleus relative to the momentum transfer direction
    std::string datafile;
    
    
};


#endif /* ipglasma_hpp */
