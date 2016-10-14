//
//  modified_ipsat.cpp
//  SubNucleon diffraction
//
//  Created by Heikki Mantysaari on 8/16/16.
//  Copyright © 2016 Heikki. All rights reserved.

#include <iostream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

// IPsat 2012
extern "C" {
    double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};



using namespace std;

// Evaluate ipsat at given x,r, anomalous dimension,  and with given proton density tp
double IpsatWithAnomalousDimension(double xbj, double r, double gamma, double tp)
{
    // read amplitude at b=0
    int IPSAT12_PAR = 1;    // m_c=1.4 GeV
    double tmpb=0;
    double n = dipole_amplitude_(&xbj, &r, &tmpb, &IPSAT12_PAR)/2.0;
    if (n >= 1)
    {
        //cout << "Dipole " << r <<", amplitude=1" << endl;
        return 1.0;
    }
    
    
    //cout << xbj << " " << r << " " << n << " " << 2.0*n - n*n << endl;
    // get exponent
    double exponent = -log(1-n);
    
    // remove T(b=0)
    double bp=4;
    exponent /= exp(0)/(2.0*M_PI*bp);
    
    // add desired density
    exponent *= tp;
    
    // effect of anomalous dimension
    exponent /= r*r;
    exponent *= pow(r, 2.0*gamma);
    
    // to adjoint rep
    n = 1.0 - exp(-exponent);
    //return n;
    
    double adjn = 2.0*n - n*n;
    //cout << "modified " << n << " " << adjn << endl;
    return adjn;
    
}

struct roothelper
{
    double tp;
    double xp;
    double gamma;
};

double roothelperf(double r, void* p)
{
    roothelper *par = (roothelper*)p;
    double nadj = IpsatWithAnomalousDimension(par->xp, r, par->gamma, par->tp);
    //cout << r << " " << nadj - (1.0 - exp(-0.5)) << endl;
    return nadj - (1.0 - exp(-0.5));
}

int main()
{
    
    const double deltaTp = 1.09648;
    
    for (int tind=0; tind<100; tind++)
    {
        for (double y=0; y<=10.75; y+=0.25)
        {
            double t = 0.0001*pow(deltaTp, tind);
            
            gsl_function F;
            F.function = &roothelperf;
            roothelper par;
            par.xp = 0.01*exp(-y);
            par.tp=t;
            par.gamma=1.2;
            F.params=&par;
            
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;
            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            double x_lo=1e-5; double x_hi = 1000;
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            int iter=0;
            int status;
            ; int max_iter = 100;
            double r;
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                r = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi,
                                                 0, 0.001);
                
            }
            while (status == GSL_CONTINUE && iter < max_iter);
            gsl_root_fsolver_free(s);
            
            cout << y << " " << t << " " << 2.0/(r*r) << endl;
            
        }
        
    }
    
    
    
}