#include <iostream>
#include <complex>
#include "vector.hpp"
#include <gsl/gsl_rng.h>
#include "gauss_boost.hpp"
#include "ipsat_proton.hpp"
#include "diffraction.hpp"
#include "subnucleon_config.hpp"
//#include "src/wilsonline.hpp"

using namespace std;
/*
int main()
{
std::vector < std::vector < std::complex<double> > > data;

vector <complex<double> > l1; 
vector <complex< double > > l2;
vector <complex < double > > l3;

l1.push_back(1);
l1.push_back(2 + 1i);
l1.push_back(-1-4i);

}
*/

#include "unit_test_framework.hpp"

gsl_rng* global_rng;

// TEST takes in one argument: the name of the test case.
// Note that the name of the test case must be a valid function name in C++.
TEST(vector_class) {
    global_rng = gsl_rng_alloc(gsl_rng_default);
    Vec v1(1,2);
    Vec v2(5,-1);
    ASSERT_ALMOST_EQUAL(v1*v2, 3,1e-7);  // If spam != 42, this test case will fail
}

TEST(forward_jpsi_amplitude_ipsat_mzwf)
{
    global_rng = gsl_rng_alloc(gsl_rng_default);
    BoostedGauss wf("gauss-boosted_mzsat.dat");
    Ipsat_Proton proton(MZSAT);
    proton.SetProtonWidth(0);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();

    Diffraction diff(proton, wf);
    
    double xp=1e-3; double Qsqr=10; double t=0.1;
    MCINTPOINTS=1e7;
    ASSERT_ALMOST_EQUAL(diff.ScatteringAmplitude(xp, Qsqr, t, T),0.04884376518,1e-5); 

}

// DO NOT REMOVE
// Generates a main() function that runs all of your tests.
// Note: Some versions of g++ incorrectly produce a warning about empty
// statements when using the -pedantic flag. Therefore, we will not put
// a semicolon after the TEST_MAIN() macro.
TEST_MAIN()
