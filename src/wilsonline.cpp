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

using std::cerr;
using std::endl;


WilsonLine::WilsonLine()   // Initialize everything to 0
{
    // By default 3x3 zero matrix
    for (int i=0; i<3; i++)
    {
        std::vector< std::complex<double> > tmpvec;
        for (int j=0; j<3; j++)
        {
            tmpvec.push_back(0);
        }
        data.push_back(tmpvec);
    }
}
WilsonLine::WilsonLine(std::vector< std::vector< std::complex<double> > >  &d)
{
    data=d;
}

void WilsonLine::Set(int row, int column, std::complex<double> value)
{
    if (row >= data.size())
    {
        cerr << "Invalid row index " << row << " num of rows in matrix " << data.size() << endl;
        return;
    }
    if (column >= data[row].size())
    {
        cerr << "Invalid column index " << column << " num of rows in matrix " << data[row].size() << endl;
        return;
    }
    
    data[row][column]=value;

}

WilsonLine WilsonLine::operator*(WilsonLine& w)
{
    std::vector< std::vector< std::complex< double> > > newdata;
    
    
    if (Size() != w.Size())
    {
        throw "WrongMatrixSize";
    }
    
    for (int row=0; row < Size(); row++)
    {
        std::vector< std::complex< double> > tmprow;
        for (int col=0; col < Size(); col++)
        {
            // Calculate element (row,column)
            // row:th line times col:th column
            std::complex<double> tmp=0;
            for (int i=0; i<Size(); i++)
            {
                tmp = tmp + Element(row, i) * w.Element(i, col);
            }
            tmprow.push_back(tmp);
        }
        newdata.push_back(tmprow);
    }
    WilsonLine res(newdata);
    return res;
}

WilsonLine WilsonLine::operator*(std::complex<double> t)
{
    std::vector< std::vector< std::complex< double> > > newdata;
    
    
    for (int row=0; row < Size(); row++)
    {
        std::vector< std::complex< double> > tmprow;
        for (int col=0; col < Size(); col++)
        {
            // Calculate element (row,column)
            tmprow.push_back(t*Element(row, col));
        }
        newdata.push_back(tmprow);
    }
    WilsonLine res(newdata);
    return res;
}


WilsonLine WilsonLine::operator+(WilsonLine& w)
{
    std::vector< std::vector< std::complex< double> > > newdata;
    
    
    if (Size() != w.Size())
    {
        throw "WrongMatrixSize";
    }
    
    for (int row=0; row < Size(); row++)
    {
        std::vector< std::complex< double> > tmprow;
        for (int col=0; col < Size(); col++)
        {
            // Calculate element (row,column)
            tmprow.push_back(Element(row, col) + w.Element(row,col));
        }
        newdata.push_back(tmprow);
    }
    WilsonLine res(newdata);
    return res;
}




WilsonLine WilsonLine::ComplexConjugate()
{
    std::vector< std::vector< std::complex<double> > > tmpdata = data;
    
    for (int row = 0; row < tmpdata.size(); row++)
    {
        for (int col = 0; col < tmpdata[0].size(); col++)
        {
            tmpdata[row][col] = std::conj( data[row][col] );
        }
    }
    
    WilsonLine wl(tmpdata);
    return wl;
}

WilsonLine WilsonLine::Transpose()
{
    std::vector< std::vector< std::complex<double> > > tmpdata = data;
    
    for (int row = 0; row < tmpdata.size(); row++)
    {
        for (int col = 0; col < tmpdata[0].size(); col++)
        {
            tmpdata[row][col] = data[col][row];
        }
    }
    
    WilsonLine wl(tmpdata);
    return wl;
}

WilsonLine WilsonLine::HermitianConjugate()
{
    WilsonLine res = ComplexConjugate().Transpose();
    return res;
}

std::complex<double> WilsonLine::Trace()
{
    std::complex<double> res;
    for (int i=0; i<Size(); i++)
        res += Element(i,i);
    return res;
}

int WilsonLine::Size()
{
    return data.size();
}

std::complex<double> WilsonLine::Element(int row , int col)
{
    return data[row][col];
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
    data.clear();
    
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
    
    
}

void WilsonLine::InitializeAsIdentity()
{
    if (data.size() != 3)
    {
        std::cerr << "WilsonLine size is on 3x3!" << std::endl;
        return;
    }
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            if (i==j)
                data[i][j]=1.0;
            else
                data[i][j]=0;
        }
    }
}

#ifdef WILSONLINE_GSL
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
    gsl_matrix_complex* m = GetGslMatrixl();
    gsl_matrix_complex* exp = gsl_matrix_complex_alloc(3,3);
    my_gsl_complex_matrix_exponential(exp,m,3);
    WilsonLine result;
    result.InitializeAsGslMatrix(exp);

    gsl_matrix_complex_free(m);
    gsl_matrix_complex_free(exp);
    
    return result;
}


gsl_matrix_complex* WilsonLine::GetGslMatrixl()
{
    gsl_matrix_complex *m = gsl_matrix_complex_alloc(3,3);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            gsl_complex c = gsl_complex_rect(Element(i,j).real(), Element(i,j).imag());
            gsl_matrix_complex_set(m, i, j, c);
        }
    }
    return m;
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
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            gsl_complex tmp;
            tmp = gsl_matrix_complex_get(m, i, j);
            std::complex<double> c(GSL_REAL(tmp),GSL_IMAG(tmp));
            data[i][j] = c;
        }
    }
}


#endif  // gsl;

