/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#ifndef subnucleon_config_h
#define subnucleon_config_h

const int ZINT_INTERVALS = 15;
const double ZINT_RELACCURACY = 0.3;
const double MCINTACCURACY = 0.2;
extern int MCINTPOINTS ;

extern bool REAL_PART;  // Calculate real part of the amplitude

enum MCINTEGRAL
{
    VEGAS,
    MISER
};

const MCINTEGRAL MCINT = VEGAS;

enum PROCESS
{
    COHERENT,
    INCOHERENT
};

#endif /* subnucleon_config_h */
