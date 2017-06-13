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

class IPGlasma : public DipoleAmplitude {
public:
    
    IPGlasma(std::string fname);
    
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2]);
    double AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] );

    WilsonLine& GetWilsonLine( double x, double y); // Find Wilson line that corresponds to the coordinate
    
    std::string InfoStr();
    
    double MinX();
    double MaxX();
    
    double XStep(); // Grid spacing in x
    
    std::vector<double> &GetXCoordinates();

    void SetSchwinger(bool s, double rc=0);
    

private:
    
    
    std::vector< double > xcoords;
    std::vector< double > ycoords;
    std::vector< WilsonLine  > wilsonlines;

    bool schwinger; // true if we use schwinger mechanism
    double schwinger_rc; // Use Schwinger for dipoles larger than schwinger_rc
    
    std::string datafile;
    
    
};


#endif /* ipglasma_hpp */
