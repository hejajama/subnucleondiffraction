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

// g++ UPC_diff_cohenrent_cuba.cpp -lcuba -lgsl -lgslcblas -lm -fopenmp -o Diff
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

double b_thetab_integrand( double *vec, size_t dim, void* par) {
    PARAMETERS *helper = (PARAMETERS*)par;
    double bmag = vec[0];
    double theta_b=vec[1];
    std::complex<double> imag(0.,1.);
    double alpha_em = 0.0073;
    double delta = std::sqrt(helper->t);
    double Bx = helper->Bmag * cos(helper->thetaB);
    double By = helper->Bmag * sin(helper->thetaB);
    double bx = bmag * cos(theta_b);
    double by = bmag * sin(theta_b);
    double mproton = 0.938;// GeV
    
    // Create the interpolation object
    gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, lb, lthetab);
    // Initialize the interpolation object with the input data
    gsl_interp2d_init(interp, mag_b_table, theta_b_table, Arealmean, lb, lthetab);
    // Evaluate the interpolated function at a new point
    double Arealmean_interp = gsl_interp2d_eval(interp, mag_b_table, theta_b_table, Arealmean, bmag, theta_b, NULL, NULL);
    gsl_interp2d_free(interp); // Need to free the interpolation object
    
    // Create the interpolation object
    interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, lb, lthetab);
    // Initialize the interpolation object with the input data
    gsl_interp2d_init(interp, mag_b_table, theta_b_table, Aimagmean, lb, lthetab);
    // Evaluate the interpolated function at a new point
    double Aimagmean_interp = gsl_interp2d_eval(interp, mag_b_table, theta_b_table, Aimagmean, bmag, theta_b, NULL, NULL);
    gsl_interp2d_free(interp); // Need to free the interpolation object
    
    double Ftilde_part = helper->Z_nuclear * std::sqrt(alpha_em) / 2./M_PI /M_PI; 
    
    std::complex<double> niAtilde(Arealmean_interp, Aimagmean_interp);
    std::complex<double> M0_part = niAtilde * Ftilde_part;
    double b_m_B = std::sqrt((Bx-bx)*(Bx-bx) + (By-by)*(By-by)) + 1.e-30;
    double b_p_B = std::sqrt((Bx+bx)*(Bx+bx) + (By+by)*(By+by)) + 1.e-30;
    std::complex<double> Mx =  std::exp(-imag*(bmag * delta * cos(theta_b))) * bmag * (
                                  M0_part * (bx-Bx)/b_m_B * helper->x2 * mproton * gsl_sf_bessel_Knu(1.0, mproton * helper->x2 * b_m_B)
                                + M0_part * (bx+Bx)/b_p_B * helper->x1 * mproton * gsl_sf_bessel_Knu(1.0, mproton * helper->x1 * b_p_B) *
                                  std::exp(-imag*(helper->Bmag * delta * cos(helper->thetaB)))
                              );
    std::complex<double> My =  std::exp(-imag*(bmag * delta * cos(theta_b))) * bmag * (
                                  M0_part * (by-By)/b_m_B * helper->x2 * mproton * gsl_sf_bessel_Knu(1.0, mproton * helper->x2 * b_m_B)
                                + M0_part * (by+By)/b_p_B * helper->x1 * mproton * gsl_sf_bessel_Knu(1.0, mproton * helper->x1 * b_p_B) *
                                  std::exp(-imag*(helper->Bmag * delta * cos(helper->thetaB)))
                              );
    double res = 0.0;
    if (helper->M12reim == 11) res = Mx.real();
    if (helper->M12reim == 12) res = Mx.imag();
    if (helper->M12reim == 21) res = My.real();
    if (helper->M12reim == 22) res = My.imag();
    //cout <<res<<endl;
    return res;
}

/* Define the integrand function */
int integrand_B_y_P(const int *ndim, const cubareal *x, const int *ncomp, cubareal *f, void *userdata){
    // x0: BigP,  x1 BigB, x2: thetaB, x3: y1, x4: y2,
    PARAMETERS *Phelper = (PARAMETERS*)userdata;
    double Bmag = Phelper->Bmin + x[1]*(Phelper->Bmax-Phelper->Bmin);
    double thetaB = x[2]*(2.*M_PI);
    double y1 = Phelper->ymin + x[3]*(Phelper->ymax-Phelper->ymin);
    double y2 = Phelper->ymin + x[4]*(Phelper->ymax-Phelper->ymin);
    double masspion = 0.13957;// GeV
    double BW_Gamma = 0.156;//GeV, rho -> pipi GeV
    double mrho_UPC = 0.77526;// GeV
    
    double BigPmin = Phelper->Qmin*Phelper->Qmin / (2.*(cosh(std::abs(y1-y2)) + 1.)) - masspion*masspion;
    if (BigPmin < 0.2) BigPmin = 0.2; // GeV
    double BigPmax = Phelper->Qmax*Phelper->Qmax / (2.*(cosh(std::abs(y1-y2)) + 1.)) - masspion*masspion;
    double BigP = BigPmin + x[0]*(BigPmax-BigPmin);
    
    double Qsquare = 2.*std::sqrt(pow(BigP*BigP + Phelper->t/4. + masspion*masspion, 2) - pow(BigP*cos(Phelper->theta_BigP), 2) * Phelper->t) * cosh(y1-y2) + 2.*(BigP*BigP - Phelper->t/4.) + 2.*masspion*masspion;
    double Bjorkonx1 = sqrt((BigP*BigP+masspion*masspion))/Phelper->rootsnn * (exp(y1) + exp(y2));
    double Bjorkonx2 = sqrt((BigP*BigP+masspion*masspion))/Phelper->rootsnn * (exp(-1.*y1) + exp(-1.*y2));

    double Mxreal, Mximag, Myreal, Myimag;
    // Do the b and theta_b integration first 
    double *lower, *upper;
    lower = new double[2];
    upper = new double[2];
    lower[0] = lower[1] = 0.;
    upper[0] = 60.;      // Max b
    upper[1] = 2.0*M_PI; // Max theta_b
    double result, error;
    const int MCINTPOINTS = 10000;
    gsl_rng* global_rng = gsl_rng_alloc(gsl_rng_default);
    
    PARAMETERS helper11;
    gsl_monte_function F11;
    F11.f = &b_thetab_integrand;
    F11.dim = 2;
    helper11.Bmag = Bmag;
    helper11.thetaB = thetaB;
    helper11.t = Phelper->t;
    helper11.x1 = Bjorkonx1;
    helper11.x2 = Bjorkonx2;
    helper11.Z_nuclear = Phelper->Z_nuclear;
    helper11.M12reim = 11;
    F11.params = &helper11;
    gsl_monte_vegas_state *s11 = gsl_monte_vegas_alloc(F11.dim);
    gsl_monte_vegas_integrate(&F11, lower, upper, F11.dim, MCINTPOINTS/50, global_rng, s11, &result, &error);
    int iter=0;
    do {
        iter++;
        gsl_monte_vegas_integrate(&F11, lower, upper, F11.dim, MCINTPOINTS/5, global_rng, s11, &result, &error);
    } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s11) - 1.0) > 0.5);
    gsl_monte_vegas_free(s11);
    Mxreal = result;
    
    PARAMETERS helper12;
    gsl_monte_function F12;
    F12.f = &b_thetab_integrand;
    F12.dim = 2;
    helper12.Bmag = Bmag;
    helper12.thetaB = thetaB;
    helper12.t = Phelper->t;
    helper12.x1 = Bjorkonx1;
    helper12.x2 = Bjorkonx2;
    helper12.Z_nuclear = Phelper->Z_nuclear;
    helper12.M12reim = 12;
    F12.params = &helper12;
    gsl_monte_vegas_state *s12 = gsl_monte_vegas_alloc(F12.dim);
    gsl_monte_vegas_integrate(&F12, lower, upper, F12.dim, MCINTPOINTS/50, global_rng, s12, &result, &error);
    iter=0;
    do {
        iter++;
        gsl_monte_vegas_integrate(&F12, lower, upper, F12.dim, MCINTPOINTS/5, global_rng, s12, &result, &error);
    } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s12) - 1.0) > 0.5);
    gsl_monte_vegas_free(s12);
    Mximag = result;
    
    PARAMETERS helper21;
    gsl_monte_function F21;
    F21.f = &b_thetab_integrand;
    F21.dim = 2;
    helper21.Bmag = Bmag;
    helper21.thetaB = thetaB;
    helper21.t = Phelper->t;
    helper21.x1 = Bjorkonx1;
    helper21.x2 = Bjorkonx2;
    helper21.Z_nuclear = Phelper->Z_nuclear;
    helper21.M12reim = 21;
    F21.params = &helper21;
    gsl_monte_vegas_state *s21 = gsl_monte_vegas_alloc(F21.dim);
    gsl_monte_vegas_integrate(&F21, lower, upper, F21.dim, MCINTPOINTS/50, global_rng, s21, &result, &error);
    iter=0;
    do {
        iter++;
        gsl_monte_vegas_integrate(&F21, lower, upper, F21.dim, MCINTPOINTS/5, global_rng, s21, &result, &error);
    } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s21) - 1.0) > 0.5);
    gsl_monte_vegas_free(s21);
    Myreal = result;
    
    PARAMETERS helper22;
    gsl_monte_function F22;
    F22.f = &b_thetab_integrand;
    F22.dim = 2;
    helper22.Bmag = Bmag;
    helper22.thetaB = thetaB;
    helper22.t = Phelper->t;
    helper22.x1 = Bjorkonx1;
    helper22.x2 = Bjorkonx2;
    helper22.Z_nuclear = Phelper->Z_nuclear;
    helper22.M12reim = 22;
    F22.params = &helper22;
    gsl_monte_vegas_state *s22 = gsl_monte_vegas_alloc(F22.dim);
    gsl_monte_vegas_integrate(&F22, lower, upper, F22.dim, MCINTPOINTS/50, global_rng, s22, &result, &error);
    iter=0;
    do {
        iter++;
        gsl_monte_vegas_integrate(&F22, lower, upper, F22.dim, MCINTPOINTS/5, global_rng, s22, &result, &error);
    } while (iter < 2 or fabs( gsl_monte_vegas_chisq(s22) - 1.0) > 0.5);
    gsl_monte_vegas_free(s22);
    Myimag = result;
    
    delete lower;
    delete upper;
    
    // Calculate the C0
    double C0 = Bmag * (Mxreal * Mxreal + Mximag * Mximag + Myreal * Myreal + Myimag * Myimag);
    // C2
    double C2 = Bmag * (Mxreal * Mxreal + Mximag * Mximag - Myreal * Myreal - Myimag * Myimag);

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
    const int ndim = 5; /* number of dimensions */
    const int ncomp = 3; /* number of components */
    cout << "Starts " <<endl;
    const long long int mineval = 10000; /* minimum number of integrand evaluations */
    cout <<"step1" <<endl;
    const long long int nvec = 1; /* minimum number of integrand evaluations */
    const cubareal epsrel = 1e-3; /* relative error tolerance */
    const cubareal epsabs = 1e-3; /* absolute error tolerance */
    const int verbose = 0; /* verbosity level */
    const long long int maxeval = 10000; /* maximum number of integrand evaluations */
    const long long int nstart = 10000;
    const long long int nincrease = 10000;
    const long long int nbatch = 10000;
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

    // Read in the mean A 2-d table
    ifstream real("rho_Q2_0.0_realpart.txt"); // open the file
    ifstream imag("rho_Q2_0.0_imagpart.txt"); // open the file
    double theta_b_step = 2.*M_PI/(lthetab-1.);
    double b_step = 60./(lb -1.); // b step!
    for (int ilb = 0; ilb < lb; ilb++) {
        mag_b_table[ilb] = b_step * ilb * 1.;
    }
    for (int ithetab = 0; ithetab < lthetab; ithetab++) {
        theta_b_table[ithetab] = theta_b_step * ithetab * 1.;
    }
    for (int ilb = 0; ilb < lb; ilb++) {
        for (int ilthetab = 0; ilthetab < lthetab; ilthetab++) {
            real >> Arealmean[ilb + ilthetab * lb];
            imag >> Aimagmean[ilb + ilthetab * lb];
             //cout <<  Arealmean[ilb + ilthetab * lb] <<endl;
        }
    }
    real.close(); // close the file
    imag.close(); // close the file

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
    sprintf(output_filename,"dsigma_dt");
    ofstream output(output_filename);
    char C0_filename[128];
    sprintf(C0_filename,"C0_dt");
    ofstream C0(C0_filename);
    char C2_filename[128];
    sprintf(C2_filename,"C2_dt");
    ofstream C2(C2_filename);
    t = 0.0;//
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
            output << integral[0] << "  ";
            C0 << integral[1] << "  ";
            C2 << integral[2] << "  ";

        }
        output << endl;
        C0 << endl;
        C2 << endl;


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


