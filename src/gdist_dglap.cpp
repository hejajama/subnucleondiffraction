/*
 * Gluon distribution from DGLAP equations
 * 
 * Uses data from file xg.dat to calculate 
 * Gluedist(x,r) = Pi^2/(2*Nc) * Alphas(x,mu(r)^2) * xg(x,r)
 *
 * This code is written by Tuomas Lappi 
 * I have made only few minor changes
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "gdist_dglap.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;


#define MINR 0.00000001
const float FMGEV=5.067731;

REAL interptransf(REAL rsqr){return -log(1.0+1.0/rsqr);} // transformation from rsqr to variable usen in interpolation

DGLAPDist::DGLAPDist(std::string file): rind(0), xind(0) {
  Intialize(file);
}

DGLAPDist::DGLAPDist() : rind(0),xind(0) {
  Intialize("xg.dat");
}

void DGLAPDist::Intialize(std::string file)
{

  ifstream gluedata(file.c_str(),ios::in);
  char linebuffer[1024];
  bool info=true;
 if(!gluedata.good()){
	cerr << " Can't find gluedata xg.dat"<< endl;
	exit(67);
	}
  do{
    gluedata >> linebuffer;
//	cerr << linebuffer << endl;
	if(!strcmp(linebuffer,"mu0sqr")){
	  gluedata >> mu0sqr;
//	  cerr << " mu0sqr " << mu0sqr << endl;
	}
	if(!strcmp(linebuffer,"C")){
	  gluedata >> C;
//	  cerr << " C " << C << endl;
	}
	if(!strcmp(linebuffer,"Nrsqr")){
	  gluedata >> Nrsqr;
//	  cerr << " Nrsqr " << Nrsqr << endl;
	}
	if(!strcmp(linebuffer,"Nxbj")){
	  gluedata >> Nxbj;
//	  cerr << " Nxbj " << Nxbj << endl;
	}
	if(!strcmp(linebuffer,"begin")){
		gluedata >> linebuffer;
	    if(!strcmp(linebuffer,"data")){info=false;}
	}
	if(!strcmp(linebuffer,"r2")){info=false; cerr << " Should not get r2 data yet" << endl;}
	
  } while(info);
  if(!(mu0sqr > 0 || C > 0 || Nrsqr > 0 || Nxbj > 0)) {cerr << "Could not decipher gdistparams" << endl;exit(34);}

  // Allocate datatables

  xbjvals = new REAL[Nxbj];
  rsqrvals = new REAL[Nrsqr];
  gdistdata = new REAL[Nrsqr*Nxbj];

  for(int i=0;i<Nrsqr;i++){
    gluedata >> linebuffer;
	if(!strcmp(linebuffer,"#")){
	  gluedata >> linebuffer;
	  if(!strcmp(linebuffer,"r2")){
	    gluedata >> rsqrvals[i];
	    rsqrvals[i] /= FMGEV * FMGEV;
	  } else {cerr << "Could not decipher gdistparams (r2) " << endl;exit(35);}
    }
  }
  for(int i=0;i<Nxbj;i++){
    gluedata >> linebuffer;
	if(!strcmp(linebuffer,"#")){
	  gluedata >> linebuffer;
	  if(!strcmp(linebuffer,"xbj")){
	    gluedata >> xbjvals[i];
	  } else {cerr << "Could not decipher gdistparams (xbj) " << endl;exit(36);}
    }
  }
  for(int i=0;i<Nxbj;i++){
    for(int j=0;j<Nrsqr;j++){
		gluedata >> gdistat(i,j);
		if(!(gdistat(i,j)>0)){cerr << "Could not decipher gdistparams (data) " << endl;exit(37);}
	}
  }
  minxbj = xbjvals[Nxbj-1];  maxxbj = xbjvals[0];
  minrsqr = rsqrvals[0];  maxrsqr = rsqrvals[Nrsqr-1];
  log2maxxbj = log2(maxxbj); deltalog2xbj = (log2(maxxbj)-log2(minxbj))/(Nxbj-1);
  log2minrsqr = log2(minrsqr); deltalog2rsqr = (log2(maxrsqr)-log2(minrsqr))/(Nrsqr-1);

  //cerr << " xbj " << minxbj << " " << maxxbj << endl;
  //cerr << " rsqr " << minrsqr << " " << maxrsqr << endl;

#ifdef SPLINE
  interpdata = new gsl_interp *[Nxbj];
  rsqrlogs = new REAL[Nrsqr];
  for(int i=0;i<Nrsqr;i++){
	rsqrlogs[i] = interptransf(rsqrvals[i]);
  }

  REAL datatemp[Nrsqr];

  for(int i=0;i<Nxbj;i++){
	interpdata[i] = gsl_interp_alloc(gsl_interp_cspline,Nrsqr);
    for(int j=0;j<Nrsqr;j++){
		datatemp[j] = gdistat(i,j);
    }
	gsl_interp_init(interpdata[i],rsqrlogs,datatemp,Nrsqr);
  }

  intaccel = gsl_interp_accel_alloc();
#endif

}; 

DGLAPDist::~DGLAPDist(){
 if(!xbjvals){delete [] xbjvals;}
 if(!rsqrvals){delete [] rsqrvals;}
 if(!gdistdata){delete [] gdistdata;}

#ifdef SPLINE
  for(int i=0;i<Nxbj;i++){
    if(!interpdata[i]){	gsl_interp_free(interpdata[i]);}
  }
  if(!interpdata){ delete [] interpdata ;}
  if(!rsqrlogs){delete [] rsqrlogs;}
  gsl_interp_accel_free(intaccel);
#endif

}; 

REAL DGLAPDist::Gluedist(REAL xbj,REAL rrsqr){

  // This function assumes that the unit of rrsqr is fm^2, so let's convert
  // rrsqr GeV^-2 => rrsqr fm^2     - H. Mäntysaari
  rrsqr /= FMGEV*FMGEV;
  double rsqr;
  if(xbj<minxbj || xbj > maxxbj){cerr << "Exceeded xbj limits " << xbj << endl; exit(56);}
  if(rrsqr < MINR){rsqr = MINR;}else{rsqr = rrsqr;}

  bool cache_bad=true;
  //  Lookup indices 
  // if(rsqr < minrsqr){rind = 0;}
  // else if(rsqr > maxrsqr){rind = Nrsqr;}
  // else 

#ifndef SPLINE

  if(rind > 0 && rind < Nrsqr){
    if(rsqr < rsqrvals[rind] && rsqr > rsqrvals[rind-1]){ // cached value not in range
		cache_bad= false;
	}
  }
  if(cache_bad){
    if (rsqr < minrsqr){
	  rind = 0;
      cache_bad = false;
    } else if (rsqr > maxrsqr){
	  rind = Nrsqr;
      cache_bad = false;
    }
  }  
  if(cache_bad){
    // at this point we know rsqr is between the limits
	rind = int((log2(rsqr)-log2minrsqr)/deltalog2rsqr) + 1;	
    // Consider this as useless checking and eventually remove:
    while(rsqr > rsqrvals[rind]){rind++; cerr << "log r not working ++ "<< rind << endl;}
    while(rsqr < rsqrvals[rind-1]){rind--; cerr << "log r not working -- "<< rind << endl;}  
  }

  cerr << "Found r ind " << rind << " : " << rsqrvals[rind-1] << " < " << rsqr << " < " <<  rsqrvals[rind] << endl;
#else
  if (rsqr < minrsqr){
	rind = 0;
  } else if (rsqr > maxrsqr){
	rind = Nrsqr; 
  } else {rind = Nrsqr-1;}
#endif
  
  // For xbj we enforced limits already
  if(!(xbj > xbjvals[xind] && xbj < xbjvals[xind-1])){ // cached value not in range
	xind = int((log2maxxbj - log2(xbj))/deltalog2xbj) + 1;	
	// cerr << "log2max " << log2maxxbj << "log2x " << log2(xbj) << " deltalog " << deltalog2xbj  << " xind " 
	//				<< xind << " floor of " << endl;
    // Consider this as useless checking and eventually remove:
//    while(xbj < xbjvals[xind]){xind++; cerr << "log x not working ++ "<< xind << endl;}
//    while(xbj > xbjvals[xind-1]){xind--; cerr << "log x not working -- "<< xind << endl;}  
  }
 
  //cerr << "Found xbj ind " << xind << " : " << xbjvals[xind-1] << " > " << xbj << " > " <<  xbjvals[xind] << endl;


 // REAL deltax = xbjvals[xind-1] - xbjvals[xind];
 // REAL deltar = rsqrvals[rind] - rsqrvals[rind-1];

  // OK, now we have indices, linear inter/extrapolation. Simplest case
 if(rind > 0 && rind < Nrsqr){
#ifndef SPLINE
   return  ( (xbj-xbjvals[xind])* ( (rsqr - rsqrvals[rind-1]) * gdistat(xind-1,rind)
				    + (rsqrvals[rind] - rsqr) * gdistat(xind-1,rind-1) ) 
	     + (xbjvals[xind-1]-xbj)* ( (rsqr - rsqrvals[rind-1]) * gdistat(xind,rind)
					+ (rsqrvals[rind] - rsqr) * gdistat(xind,rind-1) ) 
	     )/((xbjvals[xind-1] - xbjvals[xind])*(rsqrvals[rind] - rsqrvals[rind-1]));
#else
//	cout << "xind " << xind << endl << " rsqr " << rsqr <<  "rind " << rind << endl;
//	cout << " data at rinds " << gdistat(xind,rind-1) << " " << gdistat(xind,rind) << endl;
//	cout << " spline " << gsl_interp_eval(interpdata[xind],rsqrlogs,gdistdata+Nrsqr*(xind),interptransf(rsqr),intaccel) << endl;
//    cout << " old interp " << ( (rsqr - rsqrvals[rind-1]) * gdistat(xind,rind)
//					+ (rsqrvals[rind] - rsqr) * gdistat(xind,rind-1) ) /   (rsqrvals[rind] - rsqrvals[rind-1]) << endl;

    return ( (xbj-xbjvals[xind]) * gsl_interp_eval(interpdata[xind-1],rsqrlogs,gdistdata+Nrsqr*(xind-1),interptransf(rsqr),intaccel)  
			+   (xbjvals[xind-1]-xbj)* gsl_interp_eval(interpdata[xind],rsqrlogs,gdistdata+Nrsqr*(xind),interptransf(rsqr),intaccel)  ) /
					(xbjvals[xind-1] - xbjvals[xind]);
#endif   
 } else if (rind == 0){
   // extrapolate linearly in r^2 to small r from first 4 points
//   REAL extraprsqr = rsqr;
//   REAL extraprsqr0 = minrsqr;
//   REAL extraprsqr1 = rsqrvals[1];
//   REAL extraprsqr2 = rsqrvals[2];
//   REAL extraprsqr3 = rsqrvals[3];

   // try extrapolating in  log(1/r)
   REAL extraprsqr = interptransf(rsqr);
   REAL extraprsqr0 = interptransf(minrsqr);
   REAL extraprsqr1 = interptransf(rsqrvals[1]);
   REAL extraprsqr2 = interptransf(rsqrvals[2]);
   REAL extraprsqr3 = interptransf(rsqrvals[3]);

   REAL rsum1a = (1.0/3.0)*(gdistat(xind-1,0))*
     ( (extraprsqr1 - extraprsqr)/(extraprsqr1 - extraprsqr0)
       +(extraprsqr2 - extraprsqr)/(extraprsqr2 - extraprsqr0)
       +(extraprsqr3 - extraprsqr)/(extraprsqr3 -  extraprsqr0)
       );
   REAL rsum1b = (1.0/3.0)*(extraprsqr - extraprsqr0)*
     ( (gdistat(xind-1,1))/(extraprsqr1 - extraprsqr0)
       + (gdistat(xind-1,2))/(extraprsqr2 - extraprsqr0)
       + (gdistat(xind-1,3))/(extraprsqr3 - extraprsqr0)
       );
   REAL rsum2a = (1.0/3.0)*(gdistat(xind,0))*
     ( (extraprsqr1 - extraprsqr)/(extraprsqr1 - extraprsqr0)
       +(extraprsqr2 - extraprsqr)/(extraprsqr2 - extraprsqr0)
       +(extraprsqr3 - extraprsqr)/(extraprsqr3 -  extraprsqr0)
       );
   REAL rsum2b = (1.0/3.0)*(extraprsqr - extraprsqr0)*
     ( (gdistat(xind,1))/(extraprsqr1 - extraprsqr0)
       + (gdistat(xind,2))/(extraprsqr2 - extraprsqr0)
       + (gdistat(xind,3))/(extraprsqr3 - extraprsqr0)
       );
    return  ( (xbj-xbjvals[xind])* ( rsum1a + rsum1b ) 
	     + (xbjvals[xind-1]-xbj)* ( rsum2a + rsum2b   ) 
	     )/(xbjvals[xind-1] - xbjvals[xind]);
  
 } else { // rind == Nrsqr
   // extrapolate linearly in 1/rsqr to large r from last 4 points
   // formulas same as above, except replace rsqr by inverse
   // and count from end
//   REAL invrsqr = 1.0 / rsqr;
//   REAL invrsqr0 = 1.0 / rsqrvals[Nrsqr - 1];
//   REAL invrsqr1 = 1.0 / rsqrvals[Nrsqr - 2];
//   REAL invrsqr2 = 1.0 / rsqrvals[Nrsqr - 3];
//   REAL invrsqr3 = 1.0 / rsqrvals[Nrsqr - 4];
 
   // Extrapolate linearly in ln(r) instead
   REAL invrsqr = interptransf(rsqr);
   REAL invrsqr0 = interptransf(rsqrvals[Nrsqr - 1]);
   REAL invrsqr1 = interptransf(rsqrvals[Nrsqr - 2]);
   REAL invrsqr2 = interptransf(rsqrvals[Nrsqr - 3]);
   REAL invrsqr3 = interptransf(rsqrvals[Nrsqr - 4]);


   REAL rsum1a = (1.0/3.0)*(gdistat(xind-1,Nrsqr-1))*
     ( (invrsqr1 - invrsqr)/(invrsqr1 - invrsqr0)
       +(invrsqr2 - invrsqr)/(invrsqr2 - invrsqr0)
       +(invrsqr3 - invrsqr)/(invrsqr3 -  invrsqr0)
       );
   REAL rsum1b = (1.0/3.0)*(invrsqr - invrsqr0)*
     ( (gdistat(xind-1,Nrsqr-2))/(invrsqr1 - invrsqr0)
       + (gdistat(xind-1,Nrsqr-3))/(invrsqr2 - invrsqr0)
       + (gdistat(xind-1,Nrsqr-4))/(invrsqr3 - invrsqr0)
       );
   REAL rsum2a = (1.0/3.0)*(gdistat(xind,Nrsqr-1))*
     ( (invrsqr1 - invrsqr)/(invrsqr1 - invrsqr0)
       +(invrsqr2 - invrsqr)/(invrsqr2 - invrsqr0)
       +(invrsqr3 - invrsqr)/(invrsqr3 -  invrsqr0)
       );
   REAL rsum2b = (1.0/3.0)*(invrsqr - invrsqr0)*
     ( (gdistat(xind,Nrsqr-2))/(invrsqr1 - invrsqr0)
       + (gdistat(xind,Nrsqr-3))/(invrsqr2 - invrsqr0)
       + (gdistat(xind,Nrsqr-4))/(invrsqr3 - invrsqr0)
       );
   
   return  ( (xbj-xbjvals[xind])* ( rsum1a + rsum1b ) 
	     + (xbjvals[xind-1]-xbj)* ( rsum2a + rsum2b   ) 
	     )/(xbjvals[xind-1] - xbjvals[xind]);
 }
 
 // Should not get here  
 return 0;
};


ostream& operator<<(ostream& os,const DGLAPDist& ic){
	return os << " GDist object " ;
}


