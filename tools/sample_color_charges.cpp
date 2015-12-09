/*
 * Sample color charges from saturation scale, output as Wilson lines
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "sample_color_charges.hpp"
#include "../src/ipsat_proton.hpp"
#include "../src/wilsonline.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <iostream>
gsl_rng* global_rng;

using namespace std;

int main(int argc, char* argv[])
{
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    Sampler sampler;
    WilsonLine l;
    for (int i=1; i<=8; i++)
    {
        double rho =sampler.RandomColorCharge(0,0,0.01);
        //cout << i << " " << rho << endl;
        WilsonLine ta; ta.InitializeAsGenerator(i);
        WilsonLine tmpline; tmpline = ta*rho;
        l = l + tmpline;
    }
   // cout << l << endl;
    
    sampler.FillColorCharges(0.01);
    
    return 0;
}

void Sampler::FillColorCharges(double xbj)
{
    // Fill coordinate grid
    rho.clear();
    coordinates.clear();
    double maxr = 5;
    int xpoints = 100;
    double step = (2.0*maxr / xpoints);
    for (double y = -maxr; y<= maxr; y+=step)
    {
        coordinates.push_back(y);
        std::vector< std::vector < double > > tmpvec;
        for (double x = -maxr; x<= maxr; x+= step)
        {
            std::vector<double> tmprho;
            for (int a=1; a<=8; a++)
            {
                tmprho.push_back( RandomColorCharge(x, y, xbj) );
            }
          tmpvec.push_back(tmprho);
        }
        rho.push_back(tmpvec);
    } 


    
}



/* 
 * Sample one random color charge at the given point in transverse space
 */
double Sampler::RandomColorCharge(double x, double y, double xbj)
{
    double qs = SaturationScale(x, y, xbj);
    // Q_s = K g^2 mu, which gives
    // Qs^2 = K^2 (4pi as) g^2 mu^2
    // We ened g^2mu^2, use alphas=0.2 here and N=0.6
    double as = 0.2; double K = 0.6;
    double gmusqr = qs*qs / ( K*K*4.0*M_PI*as);
    
    // Sample color charge from the ColorChargeDistribution
    double rho;
    do{
        rho = gsl_rng_uniform(global_rng) ;
    } while (gsl_rng_uniform(global_rng) > ColorChargeDistribution(rho, gmusqr));

    return rho;
}

Sampler::Sampler()
{
    proton.SetProtonWidth(0.000001);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();
    Ny=100;


}


/*
 * Color charge distribution is Gaussian and local in color and transverse coordinate
 * Width is obtained from saturation scale
 * As a width this uses width/Ny, and widht = g^2 mu^2
 */
double Sampler::ColorChargeDistribution(double rho, double width)
{
    width = width / Ny;
    return std::exp(-rho*rho / (2.0 * width));    // Ok?
    
}


/*
 * Solve saturation scale defined as N(r^2=2/Qs^2) = 1 - exp(-1/2)
 */

struct SatscaleHelper
{
    double xbj;
    double x,y;
    Ipsat_Proton* proton;
};
double SatscaleHelperf(double r, void* p)
{
    SatscaleHelper *helper = (SatscaleHelper*) p;
    Vec q1(helper->x-0.5*r, helper->y);
    Vec q2(helper->x+0.5*r, helper->y);
    // Now |q1-q2|=r, 0.5(q1+q2)=b
    double res =helper->proton->Amplitude(helper->xbj, q1, q2) - (1.0 - std::exp(-0.5));
    return res;
    
}

double Sampler::SaturationScale(double x, double y, double xbj)
{
    const int MAXITER = 100;
    const double ACC = 0.0001;
    gsl_function f;
    f.function = &SatscaleHelperf;
    SatscaleHelper par;
    par.xbj = xbj;
    par.x = x;
    par.y=y;
    par.proton = &proton;
    f.params = &par;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &f, 1e-5, 100);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ACC);
    } while (iter <= MAXITER and status == GSL_CONTINUE);
    
    double satr = gsl_root_fsolver_root(s);
    
    double qs = std::sqrt( 2.0) / satr;
    
    gsl_root_fsolver_free(s);
    return qs;
    
}