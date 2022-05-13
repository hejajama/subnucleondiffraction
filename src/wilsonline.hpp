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


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

class WilsonLine
{
public:
    WilsonLine();   // Initialize everything to 0
    WilsonLine(std::vector< std::vector< std::complex<double> > >  &d);
    WilsonLine(gsl_matrix_complex *m);
    WilsonLine(const WilsonLine &m);
    
    ~WilsonLine();

    WilsonLine operator*(WilsonLine& w);
    WilsonLine operator*(std::complex<double> t);
    WilsonLine operator+(WilsonLine& w);
    WilsonLine& operator=(const WilsonLine& w);
    

    
    // Multiplies this by w^\dagger, returns the product
    // Fast in BLAS
    WilsonLine MultiplyByHermitianConjugate(WilsonLine& w);
    
    
    WilsonLine ComplexConjugate();
    WilsonLine Transpose();
    WilsonLine HermitianConjugate();
    
    void Set(int row, int column, std::complex<double> value);
    
    int Size(); // Size of NxN matrix
    std::complex<double> Element(int row, int col) const;

    std::complex<double> Trace();
    
    void InitializeAsGenerator(int a);  // Initialize as Color matrix t^a
    void InitializeAsIdentity();
    
    WilsonLine Exp();   // Calculate exponential
    
    gsl_matrix_complex* GetGslMatrix();
    
    void InitializeAsGslMatrix(gsl_matrix_complex* m);
    
private:
    gsl_matrix_complex *wline_gsl_matrix;
    
    static const int size = 3;
};


std::ostream& operator<<(std::ostream& os, WilsonLine& wl);


// Matrix exponential, calculates as writing the matrix as 2Nx2N real matrix.
// Downloaded from
void my_gsl_complex_matrix_exponential(gsl_matrix_complex *eA, gsl_matrix_complex *A, int dimx);





#endif /* wilsonline_hpp */
