/*
 * Diffraction at sub-nucleon scale
 * Wilson line / SU(3) matrix handling
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 * Standalone class, no external dependences, please!
 */

#include "wilsonline.hpp"
#include <vector>
#include <complex>
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using std::cerr;
using std::cout;
using std::endl;


WilsonLine::WilsonLine()   // Initialize everything to 0
{
    wline_gsl_matrix = gsl_matrix_complex_alloc(size,size);
    
    // By default 3x3 zero matrix
    gsl_matrix_complex_set_zero(wline_gsl_matrix);
    

}
WilsonLine::WilsonLine(std::vector< std::vector< std::complex<double> > >  &d)
{
    wline_gsl_matrix = gsl_matrix_complex_alloc(size,size);
    
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            gsl_matrix_complex_set(wline_gsl_matrix,i,j,gsl_complex_rect(d[i][j].real(),d[i][j].imag()));
        }
    }
}

WilsonLine::WilsonLine(gsl_matrix_complex* m)
{
    wline_gsl_matrix = m;
}

WilsonLine::WilsonLine(const WilsonLine &m)
{
    // Copy constructor: allocate memory for this, and copy data from m
    wline_gsl_matrix = gsl_matrix_complex_alloc(size,size);
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            std::complex<double> tmp = m.Element(i,j);
            gsl_matrix_complex_set(wline_gsl_matrix, i,j,gsl_complex_rect(tmp.real(), tmp.imag()));
        }
    }
}

WilsonLine& WilsonLine::operator=(const WilsonLine &m)
{
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            std::complex<double> tmp = m.Element(i,j);
            gsl_matrix_complex_set(wline_gsl_matrix, i,j,gsl_complex_rect(tmp.real(), tmp.imag()));
        }
    }
    return *this;
}

WilsonLine::~WilsonLine()
{
    gsl_matrix_complex_free(wline_gsl_matrix);
}

void WilsonLine::Set(int row, int column, std::complex<double> value)
{
    if (row >= size)
    {
        cerr << "Invalid row index " << row << " num of rows in matrix " << size << endl;
        return;
    }
    if (column >= size )
    {
        cerr << "Invalid column index " << column << " num of rows in matrix " << size << endl;
        return;
    }
    
    gsl_matrix_complex_set(wline_gsl_matrix,row,column,
                           gsl_complex_rect(value.real(),value.imag()));
    

}

WilsonLine WilsonLine::operator*(WilsonLine& w)
{
    gsl_matrix_complex *result = gsl_matrix_complex_alloc(size,size);
    gsl_matrix_complex_set_zero(result);
    
    
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), wline_gsl_matrix, w.GetGslMatrix(),  gsl_complex_rect(0,0), result);

    return WilsonLine(result);
}

// Multiplies this by w^\dagger, returns the product
// Fast in BLAS
WilsonLine WilsonLine::MultiplyByHermitianConjugate(WilsonLine& w)
{
    gsl_matrix_complex *result = gsl_matrix_complex_alloc(size,size);
    
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0), wline_gsl_matrix, w.GetGslMatrix(), gsl_complex_rect(0,0), result);
    
    return result;
    
}

WilsonLine WilsonLine::operator*(std::complex<double> t)
{
    std::cerr << "Multiplication by a number not implemented" << endl;
    exit(1);
}


WilsonLine WilsonLine::operator+(WilsonLine& w)
{

    std::cerr << "WilsonLine::operator+ not yet implemented for GSL matrix" << endl;
    exit(1);

    

}




WilsonLine WilsonLine::ComplexConjugate()
{

    std::cerr << "WilsonLine:ComplexConjugate not yet implemented for GSL matrix" << endl;
    exit(1);

}

WilsonLine WilsonLine::Transpose()
{

    std::cerr << "WilsonLine:Transpose not yet implemented for GSL matrix" << endl;
    exit(1);

}

WilsonLine WilsonLine::HermitianConjugate()
{
    
    WilsonLine res = ComplexConjugate().Transpose();
    return res;
}

std::complex<double> WilsonLine::Trace()
{
    gsl_complex sum;
    GSL_SET_COMPLEX(&sum, 0,0);
    
    for (unsigned int i=0; i< size; i++)
    {
        gsl_complex c = gsl_matrix_complex_get(wline_gsl_matrix, i, i);
        sum = gsl_complex_add(sum,c);
    }
    
    return std::complex<double>(GSL_REAL(sum), GSL_IMAG(sum));
}

int WilsonLine::Size()
{
    return size;
}

std::complex<double> WilsonLine::Element(int row , int col) const
{
    gsl_complex c = gsl_matrix_complex_get(wline_gsl_matrix, row, col);
    std::complex<double> cc(GSL_REAL(c), GSL_IMAG(c));
    return cc;

}
std::ostream& operator<<(std::ostream& os, WilsonLine& wl)
{
    for (int r=0; r<wl.Size(); r++)
    {
        for (int i=0; i<wl.Size(); i++)
            os << wl.Element(r,i) ;
        os << endl;
    }
    return os;
}

/* Initialize as SU(3) generator
 * a=1...8
 */
void WilsonLine::InitializeAsGenerator(int a)
{
    if (a < 1 or a >8)
    {
        std::cerr << "SU(3) generator id must be within 1...8, asked " << a << std::endl;
        return;
    }
    
    /*
    std::complex<double> imag(0,1.0);
    std::vector< std::complex<double> > row;
    switch(a)
    {
        case 1:
            row.clear();
            row.push_back(0); row.push_back(0.5); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 2:
            row.clear();
            row.push_back(0); row.push_back(-0.5*imag); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5*imag); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 3:
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(-0.5); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 4:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0.5);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 5:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-0.5*imag);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5*imag); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
        case 6:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0.5);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0.5); row.push_back(0);
            data.push_back(row);
            break;
        case 7:
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-0.5*imag);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0.5*imag); row.push_back(0);
            data.push_back(row);
            break;
        case 8:
            row.clear();
            row.push_back(1.0/(2.0*std::sqrt(3))); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(1.0/(2.0*std::sqrt(3))); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(-2.0/(2.0*std::sqrt(3)));
            data.push_back(row);
            break;
        default:
            std::cerr << "Incorrect a!" << std::endl;
    };
    */
    cerr << "Generator to gsl_matrix not yet implemented! " << endl;
    exit(1);
}

void WilsonLine::InitializeAsIdentity()
{
    
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            if (i==j)
                gsl_matrix_complex_set(wline_gsl_matrix,i,j,gsl_complex_rect(1,0));
            else
                gsl_matrix_complex_set(wline_gsl_matrix,i,j,gsl_complex_rect(0,0));
        }
    }
}


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/*
 * Calculate exp of matrix using GSL
 * Way too much memory allocation, so one should optimize if doing some
 * heavy numerics
 */
WilsonLine WilsonLine::Exp()
{
	std::cerr << "Test GSL matrix exponential before using this! WilsonLine WilsonLine::Exp() "  << std::endl;
    gsl_matrix_complex* m = GetGslMatrix();
    gsl_matrix_complex* exp = gsl_matrix_complex_alloc(3,3);
    my_gsl_complex_matrix_exponential(exp,m,3);
    WilsonLine result;
    result.InitializeAsGslMatrix(exp);

    gsl_matrix_complex_free(m);
    gsl_matrix_complex_free(exp);
    
    return result;
}


gsl_matrix_complex* WilsonLine::GetGslMatrix()
{
    return wline_gsl_matrix;
}



void my_gsl_complex_matrix_exponential(gsl_matrix_complex *eA, gsl_matrix_complex *A, int dimx)
{
    int j,k=0;
    gsl_complex temp;
    gsl_matrix *matreal =gsl_matrix_alloc(2*dimx,2*dimx);
    gsl_matrix *expmatreal =gsl_matrix_alloc(2*dimx,2*dimx);
    //Converting the complex matrix into real one using A=[Areal, Aimag;-Aimag,Areal]
    for (j = 0; j < dimx;j++)
        for (k = 0; k < dimx;k++)
        {
            temp=gsl_matrix_complex_get(A,j,k);
            gsl_matrix_set(matreal,j,k,GSL_REAL(temp));
            gsl_matrix_set(matreal,dimx+j,dimx+k,GSL_REAL(temp));
            gsl_matrix_set(matreal,j,dimx+k,GSL_IMAG(temp));
            gsl_matrix_set(matreal,dimx+j,k,-GSL_IMAG(temp));
        }
    
    gsl_linalg_exponential_ss(matreal,expmatreal,.01);
    
    double realp;
    double imagp;
    for (j = 0; j < dimx;j++)
        for (k = 0; k < dimx;k++)
        {
            realp=gsl_matrix_get(expmatreal,j,k);
            imagp=gsl_matrix_get(expmatreal,j,dimx+k);
            gsl_matrix_complex_set(eA,j,k,gsl_complex_rect(realp,imagp));
        }
    gsl_matrix_free(matreal);
    gsl_matrix_free(expmatreal);
}

void WilsonLine::InitializeAsGslMatrix(gsl_matrix_complex* m)
{
    gsl_matrix_complex_free(wline_gsl_matrix);
    wline_gsl_matrix = m;
    
}




