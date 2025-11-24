/*
 * Diffraction at sub-nucleon scale
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "subnucleon_config.hpp"
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <vector>

int MCINTPOINTS = 1e6;

bool CORRECTIONS = false;

bool FACTORIZE_ZINT = false;
// Note: There are probably factors of 4pi wrong from the
// wave functions if this is set to true


// GSL Error handler
int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
   
    // Errors related to convergence of integrals are handled when
    // gsl_integration functions are called, don't do anything with them here
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    if (gsl_errno == 14 or gsl_errno == 18 or gsl_errno == 11)
        return;

    // 15: underflows
    if (gsl_errno == 15 ) return;
    // 16: overflows
    if (gsl_errno == 16 ) return;


    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

/*
 * Str to double/int
 */
double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

std::vector<double> StrToList(std::string str)
{
    std::vector<double> tmp;
    std::stringstream buff(str);
    std::string cell;
    while (std::getline(buff, cell, ',')) {
        tmp.push_back(std::stod(cell));
    }
    return tmp;
}

