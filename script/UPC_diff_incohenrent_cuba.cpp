#include <stdio.h>
#include <cuba.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <complex>
const int lb = 501;
const int lthetab  = 301;
double Arealmean[lb*lthetab], Aimagmean[lb*lthetab], mag_b_table[lb], theta_b_table[lb];

using namespace std;
/*
   // DO the high-dimensional integration with the cuba.
   Here for the incoherent process.
   Wenbin Zhao, wenbinzhap237@gmail.com
   
   llVegas is a function in the CUBA library that can be used for multidimensional numerical integration. Here is an explanation of the parameters:
    ndim: The number of dimensions of the integral.
    ncomp: The number of components of the integrand (e.g. 1 for a scalar function, 3 for a vector function, etc.).
    integrand: A pointer to the function that computes the integrand. This function must take a const int* argument specifying the length of the vector x, a const double* argument specifying the coordinates at which to evaluate the integrand, a void* argument for passing user data, and an int argument specifying the component of the integrand to evaluate. The function must return the value of the integrand at the given point.
    userdata: A pointer to any user data that should be passed to the integrand function.
    nvec: The number of vectors to be used in the Vegas integration algorithm. It should be set to zero for the default value.
    epsrel: The relative precision of the integration.
    epsabs: The absolute precision of the integration.
    flags: A bit field of flags controlling the behavior of the integrator. It should be set to zero for the default value.
    seed: The random number seed used by the integrator.
    mineval: The minimum number of function evaluations to perform.
    maxeval: The maximum number of function evaluations to perform.
    nstart: The number of samples to take in each iteration of the Vegas algorithm.
    nincrease: The number of additional samples to take if the relative variance of the integrand is too high.
    nbatch: The number of samples to take in each iteration of the Suave algorithm.
    gridno: The number of the grid file to use for the grid integration algorithm. It should be set to zero for the default value.
    statefile: The name of a file in which to save the state of the integrator. It should be set to NULL for no saving.
    spin: A pointer to any spin information that should be passed to the integrator.
    neval: A pointer to an array that will contain the number of function evaluations performed for each component of the integrand.
    fail: A pointer to an integer that will be set to a non-zero value if the integration fails.
    integral: An array that will contain the integral of each component of the integrand.
    error: An array that will contain the estimated error of each component of the integrand.
    prob: An array that will contain the estimated probability of each component of the integrand.
*/

//  g++ UPC_diff_incohenrent_cuba.cpp -lcuba -lgsl -lgslcblas -lm -fopenmp -o Diff_incoherent
// Parameters that need to be passed to the integrand
struct PARAMETERS {
    double BigPmin;// GeV
    double BigPmax;// GeV
    double Bmin;// GeV^-1
    double Bmax;// GeV^-1
    double ymin;// GeV
    double ymax;// GeV
    double Qupcut;// GeV
    double Qlowcut;// GeV
    double rootsnn;// GeV
    int Z_nuclear;// GeV
    double t;
    double theta_BigP;
    double thetaB;
    double x1;
    double x2;
    double Bmag;
    double Qmin;
    double Qmax;
    int M12reim;
};

/* Define the integrand function */
int integrand_B_y_P(const int *ndim, const cubareal *x, const int *ncomp, cubareal *f, void *userdata){
    // x0: BigP,  x1 BigB, x2: y1, x3: y2,
    PARAMETERS *Phelper = (PARAMETERS*)userdata;
    double Bmag = Phelper->Bmin + x[1]*(Phelper->Bmax-Phelper->Bmin);
    double y1 = Phelper->ymin + x[2]*(Phelper->ymax-Phelper->ymin);
    double y2 = Phelper->ymin + x[3]*(Phelper->ymax-Phelper->ymin);
    double masspion = 0.13957;// GeV
    double BW_Gamma = 0.156;//GeV, rho -> pipi GeV
    double mrho_UPC = 0.77526;// GeV
    double massp = 0.938; //GeV
    double alpha_em = 0.0073;
    double BigPmin = Phelper->Qmin*Phelper->Qmin / (2.*(cosh(std::abs(y1-y2)) + 1.)) - masspion*masspion;
    double BigPmax = Phelper->Qmax*Phelper->Qmax / (2.*(cosh(std::abs(y1-y2)) + 1.)) - masspion*masspion;
    double BigP = BigPmin + x[0]*(BigPmax-BigPmin);
    
    double Qsquare = 2.*std::sqrt(pow(BigP*BigP + Phelper->t/4. + masspion*masspion, 2) - pow(BigP*cos(Phelper->theta_BigP), 2) * Phelper->t) * cosh(y1-y2) + 2.*(BigP*BigP - Phelper->t/4.) + 2.*masspion*masspion;
    double Bjorkonx1 = sqrt((BigP*BigP+masspion*masspion))/Phelper->rootsnn * (exp(y1) + exp(y2));
    double Bjorkonx2 = sqrt((BigP*BigP+masspion*masspion))/Phelper->rootsnn * (exp(-1.*y1) + exp(-1.*y2));

    double nX1 = Phelper->Z_nuclear*Phelper->Z_nuclear/M_PI/M_PI * alpha_em * massp * massp * Bjorkonx1 * Bjorkonx1 *
                 gsl_sf_bessel_Knu(1.0,  massp*Bjorkonx1*Bmag) * gsl_sf_bessel_Knu(1.0,  massp*Bjorkonx1*Bmag);
    double nX2 = Phelper->Z_nuclear*Phelper->Z_nuclear/M_PI/M_PI * alpha_em * massp * massp * Bjorkonx2 * Bjorkonx2 *
                 gsl_sf_bessel_Knu(1.0,  massp*Bjorkonx2*Bmag) * gsl_sf_bessel_Knu(1.0,  massp*Bjorkonx2*Bmag);
    // Calculate the C0
    double C0 = 1./2./M_PI * Bmag * (1.-gsl_sf_bessel_J0(Bmag*std::sqrt(Phelper->t))) * (nX1 + nX2);
    // C2
    double C2 = 1./2./M_PI * Bmag * gsl_sf_bessel_Jn(2.0, Bmag*std::sqrt(Phelper->t)) * (nX1 + nX2);

    double dsigmadt = 0.004031442*149.8176 * pow(BigP,3)/(pow(Qsquare-mrho_UPC*mrho_UPC, 2) + BW_Gamma*BW_Gamma*mrho_UPC*mrho_UPC) *
                      (0.5*C0 + 0.5 * C2 * cos(2.*Phelper->theta_BigP)); // frhopipi = 12.24, 1/(2pi)^3 = 0.004031442
    double C0_inte = 0.004031442*149.8176 * pow(BigP,3)/(pow(Qsquare-mrho_UPC*mrho_UPC, 2) + BW_Gamma*BW_Gamma*mrho_UPC*mrho_UPC) * C0;
    double C2_inte = 0.004031442*149.8176 * pow(BigP,3)/(pow(Qsquare-mrho_UPC*mrho_UPC, 2) + BW_Gamma*BW_Gamma*mrho_UPC*mrho_UPC) * C2;
    f[0] = dsigmadt;
    f[1] = C0_inte;
    f[2] = C2_inte;
    return 0;
}

int main(){

  // Set the number of threads to use
    int num_threads = 20;
    omp_set_num_threads(num_threads);

    double R_Nuclear = 6.37;//fm
    const double HBARC = 0.197327053; // GeV. fm
    /* Define the integration parameters */
    const int ndim = 6; /* number of dimensions */
    const int ncomp = 3; /* number of components */
    cout << "Starts " <<endl;
    const long long int mineval = 100000; /* minimum number of integrand evaluations */
    cout <<"step1" <<endl;
    const long long int nvec = 1; /* minimum number of integrand evaluations */
    const cubareal epsrel = 1e-3; /* relative error tolerance */
    const cubareal epsabs = 1e-3; /* absolute error tolerance */
    const int verbose = 0; /* verbosity level */
    const long long int maxeval = 100000; /* maximum number of integrand evaluations */
    const long long int nstart = 10000;
    const long long int nincrease = 10000;
    const long long int nbatch = 100000;
    const int gridno = 0;
    const int flags = 0; 
    const int seed = 0; 
    /* Allocate memory for the results */
    long long int neval; /* number of integrand evaluations */
    int fail; /* status flag */
    cubareal integral[ncomp]; /* integral estimates */
    cubareal error[ncomp]; /* error estimates */
    cubareal prob[ncomp]; /* CHI^2 probabilities */
    
    PARAMETERS params;

    params.Qmin    = 0.65;// GeV
    params.Qmax    = 0.9;// GeV
    params.Bmin    = 2.0 * R_Nuclear / HBARC;// GeV^-1
    params.Bmax    = 10.0 * R_Nuclear / HBARC;// GeV^-1
    params.ymin    = -1.;// GeV
    params.ymax    = 1.;// GeV
    params.Qupcut  = 0.9;// GeV
    params.Qlowcut = 0.65;// GeV
    params.rootsnn = 200.;// GeV
    params.Z_nuclear = 79;// GeV
    
    // Do the calculation
    double mint0 = 0.0001;
    double maxt3 = 0.1;
    double maxt1 = 0.0036;
    double maxt2 = 0.05;
    double tstep = 0.0002;
    double tstep2 = 0.001;
    double tstep3 = 0.01;
    int l_thetaP = 30;
    double theta_Big_P_step = M_PI/2.0/8.;
    double t;
    //output the results to file
    char output_filename[128];
    sprintf(output_filename,"dsigma_dt_incoherent");
    ofstream output(output_filename);
    
    char output_C0[128];
    sprintf(output_C0,"C0_dt_incoherent");
    ofstream C0(output_C0);
    
    char output_C2[128];
    sprintf(output_C2,"C2_dt_incoherent");
    ofstream C2(output_C2);
    
    for (t=mint0; t<=maxt3; t+=tstep) {
        output << t << "  ";
        C0 << t << "  ";
        C2 << t << "  ";
        for (double theta_Big_P = -1.*M_PI; theta_Big_P <= -M_PI/2.0+0.001; theta_Big_P+=theta_Big_P_step) {
            params.t = t;
            params.theta_BigP = theta_Big_P;
            /* Call the integrator */
            llVegas(ndim, ncomp, integrand_B_y_P, &params, nvec, epsrel, epsabs, flags, seed, mineval, maxeval, nstart, nincrease, nbatch,
                    gridno, NULL, NULL, &neval, &fail, integral, error, prob);
            /* Print the results */
            printf("Integral estimate: %e\n", integral[0]);
            printf("Error estimate: %e\n", error[0]);
            printf("Number of evaluations: %lld\n", neval);
            printf("Status flag: %d\n", fail);
            //cout << integral[0] <<endl;
            output << integral[0] << "  ";
            C0 << integral[1] << "  ";
            C2 << integral[2] << "  ";
            
        }
        output << endl;
        C0 << endl;
        C2 << endl;
        
        if (t>maxt1)
            tstep = tstep2;
        if (t>=maxt2 )
            tstep = tstep3;
    }
    output.close();
    C0.close();
    C2.close();
    return 0;
}

