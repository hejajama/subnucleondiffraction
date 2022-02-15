/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "subnucleon_config.hpp"
#include <gsl/gsl_rng.h>

int MCINTPOINTS = 1e7;

bool REAL_PART = true;

bool CORRECTIONS = false;

EXPERIMENT KINEMATICS = LHC;
FORM_FACTOR NUCLEAR_FF = POINT_CHARGE;

MCINTEGRAL MCINT = VEGAS;

bool FACTORIZE_ZINT = false;
// Note: There are probably factors of 4pi wrong from the
// wave functions if this is set to true

bool OFF_FORWARD_PHASE = true;
