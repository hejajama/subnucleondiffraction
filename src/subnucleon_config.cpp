/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "subnucleon_config.hpp"
#include <gsl/gsl_rng.h>

int MCINTPOINTS = 1e7;

bool REAL_PART = true;

bool CORRECTIONS = false;

bool FACTORIZE_ZINT = false;
// Note: There are probably factors of 4pi wrong from the
// wave functions if this is set to true
