#include <iostream>
#include <complex>
#include "vector.hpp"
#include <gsl/gsl_rng.h>
#include "gauss_boost.hpp"
#include "ipsat_proton.hpp"
#include "diffraction.hpp"
#include "subnucleon_config.hpp"
#include "virtual_photon.hpp"

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
    double eps=1e-7;
    ASSERT_ALMOST_EQUAL(v1*v2, 3,eps); 
    ASSERT_ALMOST_EQUAL((v1+v2).GetX(), 6, eps)
    ASSERT_ALMOST_EQUAL((v1+v2*(-4)).GetY(), 2+(-1)*(-4), eps)
    ASSERT_ALMOST_EQUAL(v1.Len(), std::sqrt(1*1+2*2),eps);
}

TEST(forward_jpsi_amplitude_ipsat_mzwf)
{
    global_rng = gsl_rng_alloc(gsl_rng_default);
    FACTORIZE_ZINT=false;
    BoostedGauss wf("gauss-boosted_mzsat.dat");
    Ipsat_Proton proton(MZSAT);
    proton.SetProtonWidth(0);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();

    Diffraction diff(proton, wf);
    diff.ShowVegasIterations(false);
    
    double xp=1e-3; double Qsqr=10; double t=0.1;
    MCINTPOINTS=1e6;
    ASSERT_ALMOST_EQUAL(diff.ScatteringAmplitude(xp, Qsqr, t, T),0.04884376518,1e-4); 

}

TEST(structure_function_ipsat)
{
    // Test that we reproduce structure functions as shown in https://arxiv.org/pdf/1804.05311.pdf
    global_rng = gsl_rng_alloc(gsl_rng_default);
    VirtualPhoton photon;
    Ipsat_Proton proton(MZSAT);
    proton.SetProtonWidth(0);
    proton.SetQuarkWidth(4);
    proton.InitializeTarget();
    Diffraction f2(proton, photon);
    f2.ShowVegasIterations(false);


    FACTORIZE_ZINT=true;
    double xbj=1e-3, Qsqr=2;
        
    photon.SetQuark(LIGHT, 0.03);
    MCINTPOINTS=1e4;
       
    // Use the fact that photon-proton cross section is just diffractive amplitude at t=0
    double xs_t = f2.ScatteringAmplitude(xbj, Qsqr, 0, T);
    double xs_l = f2.ScatteringAmplitude(xbj, Qsqr, 0, L);
    double structurefun = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
        
    double mc=1.3528;
    // heavy quark contribution
    photon.SetQuark(C, mc);
    double xbj_c = xbj * (1.0 + 4.0*mc*mc / Qsqr);
    double xs_t_c = 0;
    double xs_l_c = 0;
    double fl_c = 0;
    double structurefun_c = 0;
    xs_t_c = f2.ScatteringAmplitude(xbj_c, Qsqr, 0, T);
    xs_l_c = f2.ScatteringAmplitude(xbj_c, Qsqr, 0, L);
    structurefun_c = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l_c+xs_t_c);
    
        
    // note: precisio (last argument) is absolute precision, not relative
    ASSERT_ALMOST_EQUAL(structurefun+structurefun_c, 0.530119, 0.001);


   
}

// DO NOT REMOVE
// Generates a main() function that runs all of your tests.
// Note: Some versions of g++ incorrectly produce a warning about empty
// statements when using the -pedantic flag. Therefore, we will not put
// a semicolon after the TEST_MAIN() macro.
TEST_MAIN()
