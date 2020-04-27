#ifndef VEC_H
#define VEC_H

/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 

#include <iostream>
#include <vector>

typedef double REAL;

class Vec
{
    public:
        Vec();
        Vec(REAL x_, REAL y_);
        Vec(REAL x_, REAL y_, REAL z_);
        REAL GetX();
        REAL GetY();
        REAL GetZ();
        void SetX(REAL x_);
        void SetY(REAL y_);
        void SetZ(REAL z_);
        void operator+=(Vec& v);
        void operator-=(Vec& v);
        Vec& operator=(const Vec& v);
        Vec operator+(Vec& v);
        Vec operator-(Vec& v);
        void operator*=(REAL c);
        Vec operator*(REAL c);
        void Rotate3D(Vec n, double th);
        void Rotate3D(double alpha, double beta, double gamma);
    
        void Rotate2D(double angle);  // rotate x,y coords counterclockwise
    
        double operator*(Vec &v);   // dot product
    
        REAL Len();
        REAL LenSqr();
    
    private:
        REAL x,y,z;
};

std::ostream& operator<<(std::ostream& os, Vec& ic);

// Geometry functions
// Weiszfeld's algorithm to calculate geometric median (Fermat point)
Vec GeometricMedian(std::vector<Vec> &points);


#endif //VEC_H

