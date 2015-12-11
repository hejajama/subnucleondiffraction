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

// Enable GSL dependent features
#define WILSONLINE_GSL

#ifdef WILSONLINE_GSL
    #include <gsl/gsl_matrix.h>
    #include <gsl/gsl_linalg.h>
#endif
class WilsonLine
{
public:
    WilsonLine();   // Initialize everything to 0
    WilsonLine(std::vector< std::vector< std::complex<double> > >  &d);

    WilsonLine operator*(WilsonLine& w);
    WilsonLine operator*(std::complex<double> t);
    WilsonLine operator+(WilsonLine& w);
    
    WilsonLine ComplexConjugate();
    WilsonLine Transpose();
    WilsonLine HermitianConjugate();
    
    void Set(int row, int column, std::complex<double> value);
    
    int Size(); // Size of NxN matrix
    std::complex<double> Element(int row, int col);

    std::complex<double> Trace();
    
    void InitializeAsGenerator(int a);  // Initialize as Color matrix t^a
    void InitializeAsIdentity();
    
#ifdef WILSONLINE_GSL
    WilsonLine Exp();   // Calculate exponential
    
    gsl_matrix_complex* GetGslMatrixl();
    
    void InitializeAsGslMatrix(gsl_matrix_complex* m);
#endif
    
private:
    std::vector< std::vector< std::complex<double> > > data;
};


std::ostream& operator<<(std::ostream& os, WilsonLine& wl);

// gsl matrix tools
#ifdef WILSONLINE_GSL
// Matrix exponential, calculates as writing the matrix as 2Nx2N real matrix.
// Downloaded from
void my_gsl_complex_matrix_exponential(gsl_matrix_complex *eA, gsl_matrix_complex *A, int dimx);

#endif




#endif /* wilsonline_hpp */
