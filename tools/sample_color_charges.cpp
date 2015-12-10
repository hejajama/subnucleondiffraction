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
#include "fftwpp/Array.h"
#include "fftwpp/fftw++.h"
gsl_rng* global_rng;

using namespace std;
using Array::array2;

int main(int argc, char* argv[])
{
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    Sampler sampler;
    WilsonLine l;
    /*
    for (int i=1; i<=8; i++)
    {
        double rho =sampler.RandomColorCharge(0,0,0.01);
        //cout << i << " " << rho << endl;
        WilsonLine ta; ta.InitializeAsGenerator(i);
        WilsonLine tmpline; tmpline = ta*rho;
        l = l + tmpline;
    }
   // cout << l << endl;
    */
    sampler.FillColorCharges(0.01);
    sampler.FT_rho_to_k();
    return 0;
}

void Sampler::FillColorCharges(double xbj)
{
    // Fill coordinate grid
    rho.clear();
    coordinates.clear();
    double maxr = 5;
    
    double step = (2.0*maxr / xpoints);
    for (double y = -maxr; y<= maxr; y+=step)
    {
        coordinates.push_back(y);
        std::vector< std::vector < double > > tmpvec;
        std::vector< WilsonLine > rho_ta_row;
        for (double x = -maxr; x<= maxr; x+= step)
        {
            std::vector<double> tmprho;
            WilsonLine tmp_rho_t;
            for (int a=1; a<=8; a++)
            {
                double rnd_charge =RandomColorCharge(x, y, xbj);
                tmprho.push_back( rnd_charge  );
                WilsonLine ta;
                ta.InitializeAsGenerator(a);
                WilsonLine mat; mat = ta*rnd_charge;
                tmp_rho_t = tmp_rho_t + mat;    // Calculate t_a rho_a
            }
            rho_ta_row.push_back(tmp_rho_t);
            tmpvec.push_back(tmprho);
            //cout << "x: " << x << " y: " << y << " rho_a t_a: " << endl;
            //cout << tmp_rho_t << endl;
        }
        rho.push_back(tmpvec);
        rho_t.push_back(rho_ta_row);
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
    // We need g^2mu^2, use alphas=0.2 here and N=0.6
    double K = 0.6;
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
    as=0.2;
    xpoints = 50;


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
 * Fourier transfer to color charge density momentum space using FFTW
 */
void Sampler::FT_rho_to_k()
{
    unsigned int nx = coordinates.size();
    unsigned int ny = nx;
    
    cout << "Coordinates: " << endl;
    for (int i = 0; i < coordinates.size(); i++)
    {
        cout << i << "-" << coordinates[i] << "  ";
    }
    cout << endl;

    
    // Initialize grip of Wilson lines in x space
    Aplus.clear();
    for (int yind=0; yind<ny; yind++)
    {
        std::vector<WilsonLine> tmprow;
        for (int xind=0; xind<nx; xind++)
        {
            WilsonLine t;
            tmprow.push_back(t);
        }
        Aplus.push_back(tmprow);
    }
    
    // FT each component of the rho_t matrix separately
    for (int mat_y = 0; mat_y < 1; mat_y++)
    {
        for (int mat_x=0; mat_x < 1; mat_x++)
        {
            // Now create 2d grid of values of matrix[y][x]
            
            size_t align = sizeof(Complex);
            array2<Complex> data(nx,ny,align);
            fftwpp::fft2d forward(-1, data);
            fftwpp::fft2d backward(1, data);
            // Fill array, copying data, not very effective
            for (int yind=0; yind<ny; yind++)
            {
                for (int xind=0; xind<nx; xind++)
                {
                    // Debug: use Gaussian exp(-x^2)exp(-y^2)
                    data(yind,xind) = std::exp( - coordinates[xind]*coordinates[xind])*std::exp( - coordinates[yind]*coordinates[yind]);
                    //data(yind,xind) = rho_t[yind][xind].Element(mat_y, mat_x);
                }
            }
            //cout << "Color components " << mat_y << " , " << mat_x << ":  " << endl;
            //cout << data << endl;
            
            // FFT
            //cout << "Input: " << endl;
            //cout << data << endl;
            forward.fft(data);
            //cout << "Output: " << endl;
            //cout << data << endl;
            //exit(1);
            for (int ty=0; ty < ny; ty++)
            {
                //for (int tx=0; tx<nx; tx++)
                //{
                //    std::complex<double> tm = data[ty][tx];
                    //cout << ty << " " << tx << " " << tm.real() << endl;
                
                //}
                std::complex<double> tm = data[ty][ty];
                cout << ty << " " << tm << endl;
            }
            
            
            // Calculate A^+(xt) by calculating an inverse transform of g * rho(kt)/kt^2
            double g = std::sqrt(as * 4.0*M_PI); // as = g^2/(4pi)
            for (int yind=0; yind<ny; yind++)
            {
                for (int xind=0; xind<nx; xind++)
                {
                    // FFTW has exp(2pi k.x), so in practice it calculates the Fourier transfer at
                    // momenta 2pi k. To get the physical momenta we divide the components by 2pi
                    
                    ///TODO: Note that the FFTW output array is periodic, and one has to really think
                    // what is the momenta that corresponds to a given index
                    
                    double k1 = 1.0/coordinates[yind];
                    double k2 = 1.0/coordinates[xind];
                    //k1/=2.0*M_PI; k2/=2.0*M_PI;
                    double ktsqr = 1; //k1*k1 + k2*k2 + 0.1;
                    g=1;
                    if (isinf(ktsqr))
                        data(yind,xind)=0;
                    else
                    {
                        data(yind,xind) *= g/ktsqr;
                    }
                }
            }
            
            // Ft back to k space
            backward.fftNormalized(data);
            
            // Save
            for (int yind=0; yind<ny; yind++)
            {
                for (int xind=0; xind<nx; xind++)
                {
                    Aplus[yind][xind].Set(mat_y, mat_x, data[yind][xind]);
                }
            }
            
        }
    }
    
    
    // Should be ready!
    for (int i = 60; i<=100; i+=10)
        cout << Aplus[i][i].Element(0,0) << endl;
    
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