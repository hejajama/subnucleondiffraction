/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */
#ifndef subnucleon_config_h
#define subnucleon_config_h

#include <gsl/gsl_errno.h>
#include <string>
#include <vector>

const int ZINT_INTERVALS = 20;
const double ZINT_RELACCURACY = 0.0001;
const double MCINTACCURACY = 0.00001;

const double DELTA_Y = 0.1; // delta y (y=ln 1/x) used to calculate corrections

const double JPSI_MASS = 3.0969;

const double FMGEV = 5.068;
const double HBARC = 0.197327053; // GeV fm


// Globar random number generator - avoid initializing it mulitple times
#include <gsl/gsl_rng.h>

extern gsl_rng *global_rng;

#ifndef LINEINFO
#define LINEINFO __FILE__ << ":" << __LINE__
#endif

inline double SQR(double x) { return x*x; }
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);

int StrToInt(std::string str);
double StrToReal(std::string str);
std::vector<double> StrToList(std::string str);

#endif /* subnucleon_config_h */

