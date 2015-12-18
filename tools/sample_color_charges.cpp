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
#include <fstream>
#include <string>
#include <sstream>
#include "fftwpp/Array.h"
#include "fftwpp/fftw++.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
gsl_rng* global_rng;

using namespace std;
using Array::array2;

const double FMGEV = 5.067731;

int main(int argc, char* argv[])
{
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    std::vector<Sampler> samplers;
    Ipsat_Proton proton;
    proton.SetProtonWidth(0.000001);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();
    // Initialize N samplers, that is, discretize in longitudinal direction
    int nsamplers=1;
    for (int i=0; i<nsamplers; i++)
    {
        Sampler s(nsamplers, &proton);
        samplers.push_back(s);
        samplers[samplers.size()-1].FillColorCharges(0.01);
        samplers[samplers.size()-1].CalculateAplus();
        cerr <<"# done " << i << "/" << nsamplers-1 << endl;
        
    }
    cout << "# Initialized " << nsamplers << " color charge densities and calculated A^+(xt)" << endl;
    
    // Calculate Wislon lines and print them
    int points = samplers[0].GetNumOfCoordinatePoints();
    for (int yind=0; yind < points; yind++)
    {
        for (int xind=0; xind < points; xind++)
        {
            // V = prod_k exp(-i A^+_k)
            WilsonLine v;
            v.InitializeAsIdentity();
            for (int k=0; k<nsamplers; k++)
            {
                
                WilsonLine tmp = samplers[k].GetAplus(xind,yind);
                //cout << tmp << endl;
                std::complex<double> im(0,1);
                tmp = tmp*(-im);
                WilsonLine exp = tmp.Exp();
                v = v * exp;
            }
            
            
            // Print matrix, coordinates in fm
           
            cout << samplers[0].GetCoordinate(yind)/FMGEV << " " << samplers[0].GetCoordinate(xind)/FMGEV << " ";
            for (int row=0; row<3; row++)
            {
                for (int col=0; col<3; col++)
                {
                    cout << v.Element(row,col).real() << " " << v.Element(row,col).imag() << " ";;
                }
            }
            cout << endl;
        }
    }
    
    //sampler.FillColorCharges(0.01);
    //sampler.FT_rho_to_k();
    return 0;
}

void Sampler::FillColorCharges(double xbj)
{
    // Fill coordinate grid
    //rho.clear();
    rho_t.clear();
    coordinates.clear();
    double maxr = 6.8;
    double step = (2.0*maxr / xpoints);
    for (int i=0; i<xpoints; i++)
        coordinates.push_back(i);
        //coordinates.push_back(-maxr + step*i);

    // Read from file
    
    std::ifstream f("data/RhoOne.txt");
    for (int yind=0; yind<xpoints; yind++)
    {
        //coordinates.push_back(yind);
        std::vector< WilsonLine> rho_ta_row;
        for (int xind=0; xind<xpoints; xind++)
        {
            string line;
            getline(f, line);
            stringstream ss(line);
            int tmpx,tmpy;  // x and y indeces in the file
            ss >> tmpy; ss>>tmpx;
            WilsonLine tmp_rho_t;
            for (int a=1; a<=8; a++)
            {
                WilsonLine ta;
                ta.InitializeAsGenerator(a);
                WilsonLine mat;
                double rho; ss >> rho;
                mat = ta*rho;
                tmp_rho_t = tmp_rho_t + mat;
            }
            rho_ta_row.push_back(tmp_rho_t);
            
        }
        rho_t.push_back(rho_ta_row);
    }
    
    
    // Sample randomly
    /*
    for (int yind=0; yind<xpoints; yind++)
    {
        //coordinates.push_back(y);
        std::vector< std::vector < double > > tmpvec;
        std::vector< WilsonLine > rho_ta_row;
        for (int xind=0; xind<xpoints; xind++)
        {
            std::vector<double> tmprho;
            WilsonLine tmp_rho_t;
            for (int a=1; a<=8; a++)
            {
                double rnd_charge =RandomColorCharge(coordinates[xind], coordinates[yind], xbj);
                tmprho.push_back( rnd_charge  );
                WilsonLine ta;
                ta.InitializeAsGenerator(a);
                WilsonLine mat; mat = ta*rnd_charge;
                tmp_rho_t = tmp_rho_t + mat;    // Calculate t_a rho_a
            }
            rho_ta_row.push_back(tmp_rho_t);
            tmpvec.push_back(tmprho);
            //cout << y << " " << x << " " << tmp_rho_t.Element(0, 0).real() << endl;
            //cout << "x: " << x << " y: " << y << " rho_a t_a: " << endl;
            //cout << tmp_rho_t << endl;
        }
        rho.push_back(tmpvec);
        rho_t.push_back(rho_ta_row);
    }
    */
    
    //cout << rho_t[254][260] << endl;

    
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
    /*
    double rho;
    do{
        rho = gsl_rng_uniform(global_rng) ;
    } while (gsl_rng_uniform(global_rng) > ColorChargeDistribution(rho, gmusqr));
     */
    // Sample from Gaussian
    // Do as in 1502.01331
    double width = 2.0*qs*qs/Ny;
    double rho = gsl_ran_gaussian(global_rng, width);
    return rho;
}

/*
 * Color charge distribution is Gaussian and local in color and transverse coordinate
 * Width is obtained from saturation scale
 * As a width this uses width/Ny, and widht = g^2 mu^2
 */
double Sampler::ColorChargeDistribution(double rho, double width)
{
    width = width / Ny;
    width = width*width;
    return std::exp(-rho*rho / (2.0 * width));    // Ok?
    
}


Sampler::Sampler(int ny, Ipsat_Proton* proton_)
{
    
    Ny=ny;
    as=0.2;
    xpoints = 512;
    proton = proton_;

}




/*
 * Fourier transfer to color charge density momentum space using FFTW
 * Then divide by k^2 and FT back, to get A^+ in x space
 */
void Sampler::CalculateAplus()
{
    unsigned int nx = coordinates.size();
    unsigned int ny = nx;
    /*
    cout << "#Coordinates: " << endl;
    for (int i = 0; i < coordinates.size(); i++)
    {
        cout << "#" <<   i << "-" << coordinates[i] << "  ";
    }
    cout << endl;
    exit(1);*/
    

    
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
    for (int mat_y = 0; mat_y < 3; mat_y++)
    {
        for (int mat_x=0; mat_x < 3; mat_x++)
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
                    //data(yind,xind) = std::exp( - coordinates[xind]*coordinates[xind])*std::exp( - coordinates[yind]*coordinates[yind]);
                    data(yind,xind) = rho_t[yind][xind].Element(mat_y, mat_x);
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
            
            
            // Calculate A^+(xt) by calculating an inverse transform of g * rho(kt)/kt^2
            double g = std::sqrt(as * 4.0*M_PI); // as = g^2/(4pi)
            for (int yind=0; yind<ny; yind++)
            {
                for (int xind=0; xind<nx; xind++)
                {
                    // FFTW has exp(2pi k.x), so in practice it calculates the Fourier transfer at
                    // momenta 2pi k.
                    
                    // Calculate momenta
                    double delta = coordinates[1] - coordinates[0];
                    // First half of the grid points have positive momenta, starting from zero
                    // Momenta step is 1.0/(points*delta_x)
                    double k1,k2;
                    
                    
                    if (xind <= xpoints/2)
                        k1 = 1.0 / (xpoints * delta) * xind;
                    else    // At boundary goes to negative
                        k1 = -1.0/(2.0*delta) + 1.0/(xpoints*delta) * (xind - xpoints/2);
                    if (yind <= xpoints/2)
                        k2 = 1.0 / (xpoints * delta) * yind;
                    else    // At boundary goes to negative
                        k2 = -1.0/(2.0*delta) + 1.0/(xpoints*delta) * (yind - xpoints/2);
                    
                
                    
                    
                    /// TESTING with Bjoerns file:
                    g=1;
                  
                    
                    //k1 = 1.0/xpoints * (xind - xpoints/2);
                    //k2 = 1.0/xpoints * (yind - xpoints/2);
                    //k1 = -0.5 + (double)(xind)/xpoints;
                    //k2 = -0.5 + (double)(yind)/xpoints;
                    k1 *= 2.0*M_PI; k2 *=2.0*M_PI;
                    //double ktsqr = k1*k1 + k2*k2;
                    double ktsqr = 4.0*( sin(k1/2.0)*sin(k1/2.0)+sin(k2/2.0)*sin(k2/2.0)); // lattice momentum
                    
                    //if (abs(ktsqr)<0.001)
                    //    data(yind, xind) = 0;
                    //else
                        data(yind,xind) *= g/(ktsqr + 0.0023);
                    
                    if (isnan(data(yind,xind).real()) or isinf(data(yind,xind).real()))
                    {
                        cerr << "data is " << data(yind,xind) << " at yind=" << yind << " xind=" << xind << ", ktsqr = " << ktsqr << endl;
                        exit(1);
                    }
                    
                    
                    //cout << k2 << " " << k1 << " " << data(yind,xind).real() << " " << data(yind,xind).imag() << endl;
                    
                }
                //cout << endl;
            }
            //exit(1);
            
            // Ft back to k space
            backward.fftNormalized(data);
            
            //cout << "final:" << endl;
            //cout << data << endl;
            
            // Save
            for (int yind=0; yind<ny; yind++)
            {
                for (int xind=0; xind<nx; xind++)
                {
                    if (isnan(data(yind,xind).real()) or isinf(data(yind,xind).real()))
                    {
                        cerr << "After FFT back to x space data is " << data(yind,xind) << " at yind=" << yind << " xind=" << xind << endl;
                        exit(1);
                    }
                    Aplus[yind][xind].Set(mat_y, mat_x, data(yind,xind));
                }
            }
            
            
        }
    }

    
    
    
   /* // Should be ready!
    for (int i = 0; i < xpoints; i++)
    {
        for (int j=0; j<xpoints; j++)
        {
            cout << coordinates[i] << " " << coordinates[j] << " " << Aplus[i][j].Element(0,0).real() << endl;
        }
    }*/
    //for (int i = 60; i<=100; i+=10)
    //    cout << Aplus[i][i].Element(0,0) << endl;
    
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
    par.proton = proton;
    f.params = &par;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &f, 1e-6, 1000);
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

double Sampler::GetCoordinate(int ind)
{
    return coordinates[ind];
}

WilsonLine& Sampler::GetAplus(int xind, int yind)
{
    return Aplus[yind][xind];
}

int Sampler::GetNumOfCoordinatePoints()
{
    return coordinates.size();
}