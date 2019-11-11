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


double L = 2*Amplitude::FMGEV;
const double MCINTACCURACY=0.05;

struct inthelper
{
    IPGlasma *glasma;
    bool imag;
};

// Baryon profile
double T_baryon(Vec v1, Vec v2, Vec v3)
{
    const double B=4;
//    Vec center = v1+v2;
//    center = center + v3;
//    center = center*0.333;
    Vec center(0,0,0);
    v1 = v1 - center;
    v2 = v2 - center;
    v3 = v3 - center;
    
    return 1.0/std::pow(2.0*M_PI*B, 3.0) * std::exp(-(v1.LenSqr() + v2.LenSqr() + v3.LenSqr())/(2.0*B));
}

// MC vector: (bx,by,q1x,q1y,q2x,q2y,q3x,q3y)
double inthelperf_mc(double* vec, size_t dim, void* p)
{
    inthelper *par = (inthelper*)p;
    
    Vec q1(vec[0],vec[1]);
    Vec q2(vec[2],vec[3]);
    Vec q3(vec[4],vec[5]);
    Vec b(vec[6], vec[7]);
    // b independent target -> more statistics
    double q1v[2]={q1.GetX()+b.GetX(), q1.GetY()+b.GetY()};
    double q2v[2]={q2.GetX()+b.GetX(), q2.GetY()+b.GetY()};
    double q3v[2]={q3.GetX()+b.GetX(), q3.GetY()+b.GetY()};
    
    complex<double> baryon = par->glasma->BaryonOperator(0.01, q1v,q2v,q3v);
   
    if (par->imag == false) 
        return baryon.real() * T_baryon(q1,q2,q3) / std::pow(L, 2);
    else
        return baryon.imag() * T_baryon(q1,q2,q3) / std::pow(L, 2);
}




int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    int mcintpoints = StrToReal(argv[2]); //5e7;
    L = StrToReal(argv[3])*Amplitude::FMGEV;
    cout << "# Filename: " << fname << " points " << mcintpoints  << endl;

    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
    
    IPGlasma glasma(fname, 0.0, BINARY);
    glasma.SetPeriodicBoundaryConditions(true);
	
    inthelper helper;
    helper.glasma = &glasma;
    helper.imag=false;
	

    const size_t dim=8;
    gsl_monte_function fun;
    fun.params=&helper;
    fun.f = inthelperf_mc;
    fun.dim=dim;
    double min[dim] = {-L/2.0, -L/2.0,-L/2.0, -L/2.0, -L/2.0,-L/2.0,-L/2.0,-L/2.0 };
    double max[dim] = {L/2.0, L/2.0, L/2.0,  L/2.0, L/2.0,L/2.0,L/2.0,L/2.0 };
    
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
    double result,abserr; int iter=0;
    do
    {
        iter++;
        if (iter>=10)
        {
            cerr << "Mcintegral didn't converge  "<< endl;
            //return 0;
            break;
        }
        //gsl_monte_plain_integrate
        gsl_monte_miser_integrate
            (&fun, min, max, dim, mcintpoints, global_rng, s,
                               &result, &abserr);
        //cout << "# iter " << iter << " res " << result << " +/- " << abserr <<endl;
            //if (std::abs(abserr/result)>0.2)
                  //cerr << "#r=" << r << " misermc integral failed, result " << result << " relerr " << std::abs(abserr/result) << ", again.... (iter " << iter << ")" << endl;
    } while (std::abs(abserr/result)>MCINTACCURACY);


    double imag_result,imag_abserr; iter=0;
    helper.imag=true;
    do
    {
        iter++;
        if (iter>=10)
        {
            cerr << "Mcintegral didn't converge  "<< endl;
            //return 0;
            break;
        }
        gsl_monte_miser_integrate
            (&fun, min, max, dim, mcintpoints, global_rng, s,
                               &imag_result, &imag_abserr);
    } while (std::abs(imag_abserr/imag_result)>MCINTACCURACY);



    cout << result << " " << abserr << " " << imag_result << " " << imag_abserr << endl;
    //gsl_monte_plain_free (s);
    gsl_monte_miser_free(s);
    
    
   /*
    double xp=0.01;
    for (double r=0.01; r<=10; r*=1.1)
    {
        double q1[2]; double q2[2]; double q3[2];
        q1[0]=q2[0]=q3[0]=0;
        q1[1]=-r; q2[0]=r; q3[0]=r;
        std::complex<double> baryon = glasma.BaryonOperator(xp,q1,q2,q3);
        cout << r << " " << baryon.real() << " " << baryon.imag() << endl;
    }
    */
    return 0;


}
