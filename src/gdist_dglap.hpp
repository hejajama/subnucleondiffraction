#ifndef GDist_h
#define GDist_h

/*
 * Gluon distribution from DGLAP equations
 * 
 * Uses data from file xg.dat to calculate 
 * Gluedist(x,r) = Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
 *
 * This code is written by Tuomas Lappi
 * I have made only few minor changes to make this code work together
 * with my code
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */



#include <fstream>
#include <iostream>

#define SPLINE

#ifdef SPLINE
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#endif

typedef double REAL;

class DGLAPDist  {
 public: 
  DGLAPDist();
  DGLAPDist(std::string file); 
  ~DGLAPDist();

  REAL Gluedist(REAL x,REAL rsqr);

 private:
  void Intialize(std::string file);
  REAL minxbj; // these are enforced
  REAL maxxbj;
  REAL minrsqr; // this just determines whether we interpolate or extrapolate in rsqr
  REAL maxrsqr;
  REAL log2maxxbj,deltalog2xbj,log2minrsqr,deltalog2rsqr;
  int Nxbj;
  int Nrsqr;
  REAL C;
  REAL mu0sqr;
  REAL * xbjvals;
  REAL * rsqrvals;
  REAL * gdistdata;
  REAL& gdistat(int xbjind,int rsqrind){return gdistdata[Nrsqr*xbjind + rsqrind];}
  REAL& gdistat(int xbjind,int rsqrind) const {return gdistdata[Nrsqr*xbjind + rsqrind];}

#ifdef SPLINE
  gsl_interp ** interpdata;
  REAL * rsqrlogs;
  gsl_interp_accel * intaccel;
#endif

  // make these mutable to be able to cache values
  // rsqr value being asked for is between rsqrvals[rind-1] and rsqrvals[rind],
  // if outside range then rind = 0 or rind == Nrsqr
  // xbj value being asked for is between xbjvals[xind-1] and xbjvals[xind],
  mutable int rind;
  mutable int xind;

};


std::ostream& operator<<(std::ostream& os,const DGLAPDist& ic);

#endif
