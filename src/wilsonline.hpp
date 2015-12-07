/*
 * Diffraction at sub-nucleon scale
 * Wilson line / SU(3) matrix handling
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 * Standalone class, no external dependences, please!
 */

#ifndef wilsonline_hpp
#define wilsonline_hpp

#include <iostream>
#include <complex>
#include <vector>

class WilsonLine
{
public:
    WilsonLine();   // Initialize everything to 0
    WilsonLine(std::vector< std::vector< std::complex<double> > >  &d);

    WilsonLine operator*(WilsonLine& w);
    
    WilsonLine ComplexConjugate();
    WilsonLine Transpose();
    WilsonLine HermitianConjugate();
    
    void Set(int row, int column, std::complex<double> value);
    
    int Size(); // Size of NxN matrix
    std::complex<double> Element(int row, int col);

    std::complex<double> Trace();
    
    void InitializeAsGenerator(int a);  // Initialize as Color matrix t^a
    
private:
    std::vector< std::vector< std::complex<double> > > data;
};


std::ostream& operator<<(std::ostream& os, WilsonLine& wl);

#endif /* wilsonline_hpp */
