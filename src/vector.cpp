/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "vector.hpp"
#include <cmath>
#include <cstdlib>

using namespace std;
inline double SQR(double x){ return x*x; }

// **********
// Vec Class

Vec::Vec() { x=0; y=0; z=0; }
Vec::Vec(REAL x_, REAL y_) { x=x_; y=y_; z=0; }
Vec::Vec(REAL x_, REAL y_, REAL z_) { x=x_; y=y_, z=z_;}
Vec::Vec(const Vec& v) { x=v.GetX(); y=v.GetY(); z=v.GetZ(); }

void Vec::SetX(REAL x_) { x=x_; }
void Vec::SetY(REAL y_) { y=y_; }
void Vec::SetZ(REAL z_) { z=z_; }

REAL Vec::GetX() const { return x; }
REAL Vec::GetY() const { return y; }
REAL Vec::GetZ() const { return z; }

Vec& Vec::operator+=(Vec& v)
{
    x+=v.GetX();
    y+=v.GetY();
    z+=v.GetZ();
    
    return *this;
}

Vec& Vec::operator-=(Vec& v)
{
    x-=v.GetX();
    y-=v.GetY();
    z-=v.GetZ();
    
    return *this;
} 

Vec& Vec::operator=(const Vec& v)
{
    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
} 

Vec  Vec::operator+(const Vec& v)
{
    Vec sum;
    sum.SetX(x+v.GetX());
    sum.SetY(y+v.GetY());
    sum.SetZ(z+v.GetZ());
    return sum; 
}

Vec Vec::operator-(const Vec& v)
{
    Vec sum;
    sum.SetX(x-v.GetX());
    sum.SetY(y-v.GetY());
    sum.SetZ(z-v.GetZ());
    return sum; 
}

Vec& Vec::operator*=(REAL c)
{
    x*=c; y*=c; z*=c;
    
    return *this;
}

double Vec::operator*(Vec& v)
{
    return x*v.GetX() + y*v.GetY() + z*v.GetZ();
}

Vec Vec::operator*(REAL c)
{
    Vec tmp(x*c,y*c,z*c);
    return tmp;
}

REAL Vec::LenSqr()
{
    return SQR(x)+SQR(y)+SQR(z);
}

REAL Vec::Len()
{
    return sqrt(LenSqr());
}

std::ostream& operator<<(std::ostream& os, Vec& ic)
{
    return os << "Vector (" << ic.GetX() << ", " << ic.GetY() <<
        ", " << ic.GetZ() << "), |vec| = " << ic.Len() << " ";
}

void Vec::Rotate2D(double angle)
{
    // Rotate counterclokwise by angle [rad]
    double s = std::sin(angle);
    double c = std::cos(angle);
    
    double oldx = x;
    double oldy = y;
    
    x = c*oldx + s*oldy;
    y = -s*oldx + c*oldy;
    
}

/////////////////////
// Geometry
// Weiszfeld's algorithm to calculate geometric median (Fermat point)
// Ref. https://en.wikipedia.org/wiki/Geometric_median
Vec GeometricMedian(std::vector<Vec> &points)
{
    // Convergence parameters: components can change relatively/absolutely less than
    // given values
    const double ITERACCURACY_REL = 0.001;
    const double ITERACCURACY_ABS = 1e-6;
    
    const int MAXITER = 100;
    Vec y(0,0,0); // Original quess
    bool converged;
    for (unsigned int i=0; i<MAXITER; i++)
    {
        // y_{i+1} = \sum_j x_j / ||x_j - y_i||  / \sum_j 1/||x_j - y_i ||
        double normalization = 0;
        Vec newvec(0,0,0);
        for (unsigned int j=0; j<points.size(); j++)
        {
            Vec dist = points[j] - y;
            normalization += 1.0 / dist.Len();
            
            Vec tmpvec = points[j];
            tmpvec *= 1.0 / dist.Len();
            newvec += tmpvec;
        }
        
        newvec *= 1.0 / normalization;
        
        // Chec convergence
        Vec delta= newvec - y;
        converged=true;
        if ( std::abs(delta.GetX())/y.GetX() > ITERACCURACY_REL and std::abs(delta.GetX()) > ITERACCURACY_ABS )
            converged = false;
        if ( std::abs(delta.GetY())/y.GetY() > ITERACCURACY_REL and std::abs(delta.GetY()) > ITERACCURACY_ABS )
            converged = false;
        if ( std::abs(delta.GetZ())/y.GetZ() > ITERACCURACY_REL and std::abs(delta.GetZ()) > ITERACCURACY_ABS )
            converged = false;
        y = newvec;
        if (converged)
            break;
    }
    
    if (!converged)
    {
        cerr << "GeometricMedian() didn't converge!" << endl;
        cerr << "Problematic points:" << endl;
        for (unsigned int k=0; k<points.size(); k++)
            cout << points[k] << endl;
        cerr << "Best estimate: " << endl;
        cerr << y << endl;
        exit(1);
    }
    
    return y;
}


