/*
 * Wilson line / SU(3) matrix handling
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2015-2025
 */

#ifndef wilsonline_hpp
#define wilsonline_hpp

#include <iostream>
#include <complex>
#include <vector>



#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <string>

typedef unsigned int uint;

using std::cout;
using std::cerr;
using std::endl;

// Potential for slight performance boost if disable out of bound check
//#define DISABLE_OUT_OF_BOUND_CHECK

class WilsonLine 
{
public:

    WilsonLine(const std::vector<std::vector<std::complex<double>>>& mat) {
        for (int i = 0; i < NC; ++i) {
            for (int j = 0; j < NC; ++j) {
                data[i][j] = mat[i][j];
            }
        }
    }

    WilsonLine() {
        for (int i = 0; i < NC; ++i) {
            for (int j = 0; j < NC; ++j) {
                data[i][j] = 0;
            }
        }
    }
    // Note: operator() does NOT return a reference, so this would NOT work:
    // wilsonline(i,j)=1;
    // This avoids potential situation where multiple matrices might have entries
    // pointing into same memory location, but the drawback is that
    // wilsonline(i,j)=1 is perfectly valid C++, but does not do anything...
    // Defined here (as is Set) so that it can be inlined in the dipole
    // evaluations, which call it millions of times
    std::complex<double> operator()(int row, int col) const
    {
#ifndef DISABLE_OUT_OF_BOUND_CHECK
        CheckIndices(row, col);
#endif
        return data[row][col];
    }

    std::complex<double> Trace() const;

    WilsonLine operator*(WilsonLine& w);
    WilsonLine operator*(std::complex<double> t);
    WilsonLine operator+(WilsonLine& w);
    

    
    // Multiplies this by w^\dagger, returns the product
    WilsonLine MultiplyByHermitianConjugate(const WilsonLine& other) const;
    
    
    WilsonLine ComplexConjugate();
    WilsonLine Transpose();
    WilsonLine HermitianConjugate();
    
    void Set(int row, int column, std::complex<double> value)
    {
#ifndef DISABLE_OUT_OF_BOUND_CHECK
        CheckIndices(row, column);
#endif
        data[row][column] = value;
    }

    int Size(); // Size of NxN matrix
    std::complex<double> Element(int row, int col) const;
    
    void InitializeAsGenerator(int a);  // Initialize as Color matrix t^a
    void InitializeAsIdentity();
    
    WilsonLine Exp();   // Calculate exponential
private:
    static const int NC = 3;

    std::complex<double> data[NC][NC];

    // Throws std::out_of_range if (row, column) is not a valid matrix index
    static void CheckIndices(int row, int column)
    {
        if (row >= NC or row < 0)
            throw std::out_of_range("Invalid row index " + std::to_string(row) +  ", num of rows in matrix:  "  + std::to_string(NC));
        if (column >= NC or column < 0)
            throw std::out_of_range("Invalid column index " + std::to_string(column) +  ", num of columns in matrix: "  + std::to_string(NC));
    }

    

    
};


std::ostream& operator<<(std::ostream& os, WilsonLine& wl);





#endif /* wilsonline_hpp */
