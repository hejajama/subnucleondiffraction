/*
 * Sample color charges from saturation scale, output as Wilson lines
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef sample_color_charges_hpp
#define sample_color_charges_hpp
#include <vector>

#include "../src/ipsat_proton.hpp"
#include "../src/wilsonline.hpp"

/*
 * This class takes care of sampling color charges at one longitudinal coordinate index
 * The final output is the grid of matrices A^+ for each transverse coordinate
 * In order to finally get the Wilson line one should sample N_y matrices and calculate the 
 * product of exonentials
 */
class Sampler
{
public:
    Sampler(int ny, Ipsat_Proton* proton_);
    
    double SaturationScale(double x, double y, double xbj); // Solve saturation sacle at coordinate (x,y) at given xbj
    double ColorChargeDistribution(double rho, double width);
    double RandomColorCharge(double x, double y, double xbj);
    
    void FillColorCharges(double xbj); // Fill grid
    
    void CalculateAplus();
    
    WilsonLine& GetAplus(int xind, int yind);
    double GetCoordinate(int ind);
    
    int GetNumOfCoordinatePoints();
    

private:
    Ipsat_Proton *proton;
    int Ny;
    double as;  // alpha_s
    
    std::vector< std::vector< WilsonLine > > rho_t; // rho^a t^a for each coordinate
    std::vector< std::vector< WilsonLine > > Aplus; // A^+(xt) at each transverse coordinate
    
    // Grid
    std::vector<double> coordinates;    // maps index -> coordinate in GeV^-1
    std::vector< std::vector< std::vector<double> > > rho;   // Grid of 8 component color vectors
    
    int xpoints; // Number of grid points

};


#endif /* sample_color_charges_hpp */
