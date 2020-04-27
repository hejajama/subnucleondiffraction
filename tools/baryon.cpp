/* Calculate avarage dipole amplitude averaged over b
*/

#include "../src/ipglasma.hpp"
#include "../src/ipsat_proton.hpp"
#include <tools/tools.hpp>
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "../src/vector.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <tools/tools.hpp>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_rng.h>

gsl_rng* global_rng;
using namespace std;


double L = 5*Amplitude::FMGEV;
const double MAXDIST = 5 *Amplitude::FMGEV; // Max q distance from the center
const double MCINTACCURACY=0.1;
double probeB=4.0;

bool totxs = false;
bool normalize_area = true;
bool dipole = false; // compute dipole instead of baryon
bool slopeb2=false;

struct inthelper
{
    IPGlasma *glasma;
    bool imaginary_part;
};

// Baryon profile
double T_baryon(Vec v1, Vec v2)
{
    const double B=probeB;
    // At this stage v1,v2 are vectors pointing from b to quark 1 or 2
    double x = 2.0/3.0; 
    return 3.0/std::pow(2.0*M_PI*B/x, 2) * std::exp(-2*x*(v1.LenSqr() + v2.LenSqr() + v1*v2)/(2.0*B));
}
//
// MC vector: (b1x,b1y,q1x,q1y,q2x,q2y)
// Location for vec q3 is fixed by q1+q2+q3=0
// q1,q2 are measured w.r.t. impact parameter
double inthelperf_mc(double* vec, size_t dim, void* p)
{
    inthelper *par = (inthelper*)p;
    
    Vec b(vec[0],vec[1]);
    Vec q1(vec[2],vec[3]);
    Vec q2(vec[4],vec[5]);
   
    double density = T_baryon(q1,q2);
    
    double q1v[2]={q1.GetX()+b.GetX(), q1.GetY()+b.GetY()};
    double q2v[2]={q2.GetX()+b.GetX(), q2.GetY()+b.GetY()};
    double q3v[2]={b.GetX() - q1.GetX() - q2.GetX(), b.GetY()-q1.GetY()-q2.GetY()};
    
    complex<double> baryon;
    if (dipole == false)
    { 
        baryon = par->glasma->BaryonOperator(0.01, q1v,q2v,q3v);
        if (totxs)
            baryon = 1.0 + baryon/6.0;
        if (slopeb2)
            baryon *= b.LenSqr()/2.0;
    }
    else
    {
        baryon = par->glasma->Amplitude(0.01, q1v, q2v);
    }
    double normalization = std::pow(L, 2);
    if (normalize_area == false)    
        normalization=1;

    if (par->imaginary_part == false)
        return baryon.real() * density / normalization;
    else
        return baryon.imag() * density /normalization;
}



double inthelperf_dipole_mc(double* vec, size_t dim, void* p)
{
    inthelper *par = (inthelper*)p;
    
    Vec b(vec[0],vec[1]);
    Vec q1(vec[2],vec[3]);
   
    double density = 1.0/(M_PI*probeB) * std::exp(-2.*q1.LenSqr() / (2.0*probeB) );
    
    double q1v[2]={q1.GetX()+b.GetX(), q1.GetY()+b.GetY()};
    double q2v[2]={b.GetX() - q1.GetX() , b.GetY()-q1.GetY()};
    
    double dipole = par->glasma->Amplitude(0.01, q1v, q2v);
    
    double normalization = std::pow(L, 2);
    if (normalize_area == false)    
        normalization=1;

    return density*dipole/normalization;

}
int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    int mcintpoints = StrToReal(argv[2]);
    L = StrToReal(argv[3])*Amplitude::FMGEV;
    string mode = argv[4];
    probeB=StrToReal(argv[5]);
    size_t dim = 4;
    if (mode == "baryon")
    {
        dipole = false;
        dim=6;
    }
    else if (mode == "elastic_slope_baryon")
    {
        dipole=false;
        dim=6;
        slopeb2=true;
    }
    else if (mode == "dipole")
    {
        dipole = true;
        dim=4;
    }
    else
    {
        cerr << "Unknown mode " << mode << endl;
        exit(1);
    }
    cout << "# Filename: " << fname << " points " << mcintpoints  << " L = " << L/Amplitude::FMGEV << " fm" << endl;
    cout <<"# Mode: " << mode << " probe size " << probeB << " GeV^(-2)" <<  endl;
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
   

    bool periodicboundary = false;
    //IPGlasma glasma(fname, 0.01, TEXT);
    //glasma.SetPeriodicBoundaryConditions(false);
	IPGlasma glasma(fname, 0.01, BINARY);
    glasma.SetPeriodicBoundaryConditions(periodicboundary);
    totxs = true;
    normalize_area = false;

/*
    // test
    for (double th = 0; th <= 2.0*M_PI; th+=0.1)
    {
        double a = 5;
        double q1[2] = {-a,-a};
        double q2[2] = {a,-a};
        double q3[2] = {0, (std::sqrt(3)-1.0)*a};
        Vec q1v (q1[0], q1[1] ); q1v.Rotate2D(th); q1[0] = q1v.GetX(); q1[1] = q1v.GetY(); 
        Vec q2v (q2[0], q2[1] ); q2v.Rotate2D(th); q2[0] = q2v.GetX(); q2[1] = q2v.GetY();
        Vec q3v (q3[0], q3[1] ); q3v.Rotate2D(th); q3[0] = q3v.GetX(); q3[1] = q3v.GetY();

        std::complex<double> b = glasma.BaryonOperator(0.01, q1,q2,q3);
        std::complex<double> b2 = glasma.BaryonOperator(0.01, q3,q1,q2);
        std::complex<double> b3 = glasma.BaryonOperator(0.01, q1,q3,q2);
        cout << th << " " << b.real() << " " << b.imag() << endl; // " " << b2.real() << " " << b2.imag() << " " << b3.real() << " " << b3.imag() << endl; 

    }
exit(1);
*/
    cout << "# Compute totxs: " << totxs << " normalize area: " << normalize_area << " periodic boundary conditions: " << periodicboundary << endl;

    inthelper helper;
    helper.glasma = &glasma;
    helper.imaginary_part=false;
	

    gsl_monte_function fun;
    fun.params=&helper;
    double *min;
    double *max;
    if (dipole == false)
    {
        dim=6;
        min = new double[dim];
        max = new double[dim]; 
        min[0] = -L/2; min[1] = -L/2; min[2] = -MAXDIST; min[3]=-MAXDIST; min[4]=-MAXDIST; min[5]=-MAXDIST;
        max[0] = L/2; max[1] = L/2; max[2] = MAXDIST; max[3]=MAXDIST; max[4]=MAXDIST; max[5]=MAXDIST;

        fun.f = inthelperf_mc;
    }
    else
    {
        dim=4;
        min = new double[dim];
        max =new double[dim];
        min[0]=-L/2; min[1]=-L/2;  min[2] = -MAXDIST; min[3]=-MAXDIST;
        max[0] = L/2; max[1] = L/2; max[2] = MAXDIST; max[3]=MAXDIST;
        fun.f = inthelperf_dipole_mc;
   }
     
    fun.dim=dim;
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
    double result_real=0; double result_imag=0;
    
    
    /// Real part
    helper.imaginary_part=false;
    double abserr; int iter=0;
    do
    {
        iter++;
        if (iter>=10)
        {
            cerr << "Mcintegral didn't converge  "<< endl;
            break;
            //return 0;
            break;
        }
        //gsl_monte_plain_integrate
        gsl_monte_miser_integrate
            (&fun, min, max, dim, mcintpoints, global_rng, s,
                               &result_real, &abserr);
        //cout << "# iter " << iter << " res " << result << " +/- " << abserr <<endl;
            //if (std::abs(abserr/result)>0.2)
                  //cerr << "#r=" << r << " misermc integral failed, result " << result << " relerr " << std::abs(abserr/result) << ", again.... (iter " << iter << ")" << endl;
    } while (std::abs(abserr/result_real)>MCINTACCURACY);


    double imag_abserr; iter=0;
  /* 
    //// Imaginary
    iter=0;
    result_imag=0;
    helper.imaginary_part=true;
    do
    {
        iter++;
        if (iter>=3)
        {
            cerr << "Imaginary part mcintegral didn't converge, result " << result_imag << " +/- " << abserr << endl;
            break;
        }
        //gsl_monte_plain_integrate
        gsl_monte_miser_integrate
            (&fun, min, max, dim, mcintpoints, global_rng, s,
                               &result_imag, &abserr);
        cout << "# imaginary part iter " << iter << " res " << result_imag << " +/- " << abserr <<endl;
          
    } while (std::abs(abserr/result_imag)>MCINTACCURACY);
*/
    cout << result_real  << " " << result_imag << endl;
    //gsl_monte_plain_free (s);
    delete [] min;
    delete[] max; 
    gsl_monte_miser_free(s);
    
    
    return 0;


}
