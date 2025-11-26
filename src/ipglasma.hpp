/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a IPglasma nucleus
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef ipglasma_hpp
#define ipglasma_hpp

#include <string>
#include <vector>
#include "dipole.hpp"
#include "wilsonline.hpp"

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
    std::complex<double> ComplexAmplitude(double xpom, double q1[2], double q2[2]) ;
    double Amplitude(double xpom, double q1[2], double q2[2]) ;
    double AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] ) ;

    const WilsonLine& GetWilsonLine( double x, double y) const; // Find Wilson line at the given point (x,y) [GeV^-1]
    const WilsonLine& GetWilsonLine(int i) const { return wilsonlines[i]; } 

    std::string InfoStr();
    
    double MinX();
    double MaxX();
    
    double XStep(); // Grid spacing in x

    void SetPeriodicBoundaryConditions(bool s) { periodic_boundary_conditions = s; }
    
    std::vector<double> &GetXCoordinates();

    void ApplyPeriodicBoundaryConditions(double q[2]) const; 

    std::vector<int> LatticeCoordinates(double x, double y) const;
    int WilsonLineCoordinate(int  xind, int yind) const;
    

    double X(int ix) { return xcoords[ix]; }
    double Y(int iy) { return ycoords[iy]; }

    std::vector<int> LatticeCoordinates(double x, double y);
    int WilsonLineCoordinate(int  xind, int yind);
    WilsonLine& GetWilsonLine(int i) { return wilsonlines[i]; }

private:
    
    
    std::vector< double > xcoords;
    std::vector< double > ycoords;
    std::vector< WilsonLine  > wilsonlines;
    
    bool periodic_boundary_conditions;
    
    std::string datafile;

    static const int NC=3;
    
    
};


#endif /* ipglasma_hpp */
