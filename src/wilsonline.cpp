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
    
    switch(a)
    {
        case 1:
            std::vector< std::complex<double> > row;
            row.push_back(0); row.push_back(0.5); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0.5); row.push_back(0); row.push_back(0);
            data.push_back(row);
            row.clear();
            row.push_back(0); row.push_back(0); row.push_back(0);
            data.push_back(row);
            break;
    };
    
    
}


