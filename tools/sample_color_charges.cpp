/*
 * Sample color charges from saturation scale, output as Wilson lines
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "sample_color_charges.hpp"
#include "../src/ipsat_proton.hpp"
#include "../src/wilsonline.hpp"
#include "tools/tools.hpp"
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
double QS_COLOR_CHARGE_COEF = 0.34;   // S(b) = c Q_s^2
int GRIDPOINTS = 256;
//double MAXR = 12.0*FMGEV;
double MAXR = 3.0*FMGEV;
double XPOM = 0.000959089;
const int LONGITUDINAL_STEPS = 40;



int main(int argc, char* argv[])
{
    if (argc > 1)
        GRIDPOINTS = StrToInt(argv[1]);
    if (argc > 2)
        QS_COLOR_CHARGE_COEF = StrToReal(argv[2]);
    
    cout << "# Sampling random color charges" << endl;
    cout << "# Heikki Mäntysaari <mantysaari@bnl.gov>, 2015-2016" << endl;
    cout << "# Params: S(b)=c Q_s^2 with c=" << QS_COLOR_CHARGE_COEF << ", gridpoints " << GRIDPOINTS << " grid size [-" << MAXR << ", " << MAXR << "], x=" << XPOM << ", longitudinal steps " << LONGITUDINAL_STEPS << endl;
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    
    std::vector<Sampler> samplers;
    Ipsat_Proton proton;
    proton.SetProtonWidth(0.000001);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();
    // Initialize N samplers, that is, discretize in longitudinal direction
    int nsamplers=LONGITUDINAL_STEPS;
    for (int i=0; i<nsamplers; i++)
    {
        Sampler s(nsamplers, &proton);
        samplers.push_back(s);
        samplers[samplers.size()-1].FillColorCharges(XPOM);
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
            //for (int row=0; row<3; row++)
            for (int col=0; col<3; col++)
            {
                //for (int col=0; col<3; col++)
                for (int row=0; row<3; row++)
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
    rho_t.clear();
    coordinates.clear();
    double maxr = MAXR;
    double step = (2.0*maxr / (xpoints-1));
    for (int i=0; i<xpoints; i++)
        //coordinates.push_back(i);
        coordinates.push_back(-maxr + step*i);

    // Read from file
    // Note: when calculating fourier transforms we change the color charge density
    // to the lattice units. Thus, the file should be in units [rho]=GeV
    /*
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
    */
    
    // Sample randomly
    
    for (int yind=0; yind<xpoints; yind++)
    {
        //coordinates.push_back(y);
        std::vector< WilsonLine > rho_ta_row;
        for (int xind=0; xind<xpoints; xind++)
        {
            WilsonLine tmp_rho_t;
            for (int a=1; a<=8; a++)
            {
                double rnd_charge =RandomColorCharge(coordinates[xind], coordinates[yind], xbj);
                // charge density has a dimension 1/(GeV^-2), make it dimensionless
                // charge/unit cell
                //double lattice_l = coordinates[ coordinates.size() - 1 ] - coordinates[0];
                //double gev_to_lattice = lattice_l / xpoints;
                //rnd_charge *= gev_to_lattice;
                
                WilsonLine ta;
                ta.InitializeAsGenerator(a);
                WilsonLine mat; mat = ta*rnd_charge;
                tmp_rho_t = tmp_rho_t + mat;    // Calculate t_a rho_a
                //cout << rnd_charge << " ";
            }
            //cout << endl;
            rho_ta_row.push_back(tmp_rho_t);
            //cout << yind << " " << xind << " " << tmp_rho_t.Element(0, 0).real() << endl;
            //cout << "x: " << x << " y: " << y << " rho_a t_a: " << endl;
            //cout << tmp_rho_t << endl;
        }
        //cout << endl;
        //cout << endl;
        rho_t.push_back(rho_ta_row);
    }
    
    //cout << rho_t[254][260] << endl;

    
}


 
/* 
 * Sample one random color charge at the given point in transverse space
 */
double Sampler::RandomColorCharge(double x, double y, double xbj)
{
    //x=0; y=0; //tmp
    Vec tmpcoord(x,y);
    double qs = 0.5; //proton->SaturationScale(xbj,tmpcoord);

    
    // Too small density to find Q_s, put color charge to 0
    if (qs < 0)
        return 0;
    
    
    // Sample color charge from the ColorChargeDistribution
    /*
    double rho;
    do{
        rho = gsl_rng_uniform(global_rng) ;
    } while (gsl_rng_uniform(global_rng) > ColorChargeDistribution(rho, gmusqr));
     */
    // Sample from Gaussian
    // Do as in 1502.01331
    double variance = QS_COLOR_CHARGE_COEF*qs*qs/Ny;
    // Note: variance is c*Q_s^2/N_y, variance = std deviation^2
    // For the gsl_ran the 2nd argument is std deviation
    double stdev = std::sqrt(variance);
    double rho = gsl_ran_gaussian(global_rng, stdev);
    return rho;
}

/*
 * Color charge distribution is Gaussian and local in color and transverse coordinate
 * Width is obtained from saturation scale
 * As a width this uses width/Ny, and widht = g^2 mu^2
 */
double Sampler::ColorChargeDistribution(double rho, double width)
{
    cerr << "dont use me!" << endl;
    width = width / Ny;
    width = width*width;
    return std::exp(-rho*rho / (2.0 * width));    // Ok?
    
}


Sampler::Sampler(int ny, Ipsat_Proton* proton_)
{
    
    Ny=ny;
    as=0.2;
    xpoints = GRIDPOINTS;
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
                    double delta = 1; //coordinates[1] - coordinates[0];
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
                    
                
                    // Note lattice units: one unit in k space is 2pi/L
                    // Bjoerns file: L = 24 fm
                    double lattice_l = coordinates[ coordinates.size() - 1 ] - coordinates[0];
                    double lattice_a = (coordinates[1] - coordinates[0])/GRIDPOINTS;
                    //double kstep = 2.0*M_PI / (lattice_l );
                    //double m_lattice_units = 0.2 / kstep * 2.0*M_PI / xpoints; // 0.2 GeV
                    double gev_to_lattice = lattice_l / xpoints;
                    double m_lattice_units = 0.2*gev_to_lattice;
                    /// TESTING with Bjoerns file:
                    g=1;
                  
                    //k1 = -0.5 + (double)(xind)/xpoints;
                    //k2 = -0.5 + (double)(yind)/xpoints;
                    k1 *= 2.0*M_PI; k2 *=2.0*M_PI;
                    //double ktsqr = k1*k1 + k2*k2;
                    double ktsqr = 4.0*( sin(k1/2.0)*sin(k1/2.0)+sin(k2/2.0)*sin(k2/2.0)); // lattice momentum
                    
                    
                    // move to physical units
                    //double ktsqr_physical = (k1*k1 + k2*k2) / (lattice_a*lattice_a);
                    //data(yind, xind) *= g/(ktsqr_physical + 0.2*0.2);
                    
                    data(yind,xind) *= g/(ktsqr + m_lattice_units*m_lattice_units);
                    
                    
                    
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