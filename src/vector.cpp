/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "vector.hpp"
#include <tools/config.hpp>

using namespace Amplitude;

// **********
// Vec Class

Vec::Vec() { x=0; y=0; z=0; }
Vec::Vec(REAL x_, REAL y_) { x=x_; y=y_; z=0; }
Vec::Vec(REAL x_, REAL y_, REAL z_) { x=x_; y=y_, z=z_;}
//Vec::Vec(Vec& v) { x=v.GetX(); y=v.GetY(); z=v.GetZ(); }

void Vec::SetX(REAL x_) { x=x_; }
void Vec::SetY(REAL y_) { y=y_; }
void Vec::SetZ(REAL z_) { z=z_; }

REAL Vec::GetX() { return x; }
REAL Vec::GetY() { return y; }
REAL Vec::GetZ() { return z; }

void Vec::operator+=(Vec& v)
{
    x+=v.GetX();
    y+=v.GetY();
    z+=v.GetZ();
}

void Vec::operator-=(Vec& v)
{
    x-=v.GetX();
    y-=v.GetY();
    z-=v.GetZ();
} 

Vec& Vec::operator=(const Vec& v)
{
    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
} 

Vec Vec::operator+(Vec& v)
{
    Vec sum;
    sum.SetX(x+v.GetX());
    sum.SetY(y+v.GetY());
    sum.SetZ(z+v.GetZ());
    return sum; 
}

Vec Vec::operator-(Vec& v)
{
    Vec sum;
    sum.SetX(x-v.GetX());
    sum.SetY(y-v.GetY());
    sum.SetZ(z-v.GetZ());
    return sum; 
}

void Vec::operator*=(REAL c)
{
    x*=c; y*=c; z*=c;
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


