/*
 * Calculate DIS cross section
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 */

#ifndef dis_hpp
#define dis_hpp

#include "dipole.hpp"
#include "wave_function.hpp"
#include "qcd.hpp"


class DIS
{
public:
    DIS(DipoleAmplitude* dipole);
    
    double PhotonProtonCrossSection(double qsqr, double xbj, Polarization pol);
    
    double F2(double qsqr, double xbj);
    
    DipoleAmplitude* GetDipole();
private:
    DipoleAmplitude* dipole;

    
};

#endif /* dis_hpp */
