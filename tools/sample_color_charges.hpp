/*
 * Sample color charges from saturation scale, output as Wilson lines
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef sample_color_charges_hpp
#define sample_color_charges_hpp
#include <vector>

#include "../src/ipsat_proton.hpp"

class Sampler
{
public:
    Sampler();
    
    double SaturationScale(double x, double y, double xbj); // Solve saturation sacle at coordinate (x,y) at given xbj
    double ColorChargeDistribution(double rho, double width);
    double RandomColorCharge(double x, double y, double xbj);
    
private:
    Ipsat_Proton proton;
    int Ny;
    
    
    
    
    
    // Grid
    std::vector<double> coordinates;    // maps index -> coordinate in GeV^-1
    std::vector< std::vector< std::vector<double> > > rho;   // Grid of 8 component color vectors
};


#endif /* sample_color_charges_hpp */
