/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#ifndef subnucleon_config_h
#define subnucleon_config_h


const int ZINT_INTERVALS = 15;
const double ZINT_RELACCURACY = 0.3;
const double MCINTACCURACY = 0.2;
const int MCINTPOINTS = 1e7;

enum MCINTEGRAL
{
    VEGAS,
    MISER
};

const MCINTEGRAL MCINT = MISER;

enum PROCESS
{
    COHERENT,
    INCOHERENT
};

#endif /* subnucleon_config_h */
