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

bool totxs = false;
bool normalize_area = true;
bool dipole = false; // compute dipole instead of baryon

struct inthelper
{
    IPGlasma *glasma;
    bool imaginary_part;
};

// Baryon profile
double T_baryon(Vec v1, Vec v2)
{
    const double B=4.0;
//    Vec center = v1+v2;
//    center = center + v3;
//    center = center*0.333;
    Vec center(0,0,0);
    v1 = v1 - center;
    v2 = v2 - center;
    
    return 1.0/std::pow(2.0*M_PI*B, 2) * std::exp(-(v1.LenSqr() + v2.LenSqr())/(2.0*B));
}
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




int main(int argc, char* argv[])
{
    // Arguments: ipglasma filename  b  schwinger_r
    string fname = argv[1];
    int mcintpoints = StrToReal(argv[2]);
    L = StrToReal(argv[3])*Amplitude::FMGEV;
    string mode = argv[4];
    if (mode == "baryon")
        dipole = false;
    else if (mode == "dipole")
        dipole = true;
    else
    {
        cerr << "Unknown mode " << mode << endl;
        exit(1);
    }
    cout << "# Filename: " << fname << " points " << mcintpoints  << " L = " << L/Amplitude::FMGEV << " fm" << endl;
    cout <<"# Mode: " << mode << endl;
    
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off ();
   

    bool periodicboundary = true;
    //IPGlasma glasma(fname, 0.01, TEXT);
    //glasma.SetPeriodicBoundaryConditions(true);
	IPGlasma glasma(fname, 0.01, BINARY);
    glasma.SetPeriodicBoundaryConditions(periodicboundary);
    totxs = true;
    normalize_area = true;;

    cout << "# Compute totxs: " << totxs << " normalize area: " << normalize_area << " periodic boundary conditions: " << periodicboundary << endl;

    inthelper helper;
    helper.glasma = &glasma;
    helper.imaginary_part=false;
	

    size_t dim=6;
    gsl_monte_function fun;
    fun.params=&helper;
    fun.f = inthelperf_mc;
    fun.dim=dim;
    double min[6] = {-L/2.0, -L/2.0, -MAXDIST,-MAXDIST,-MAXDIST,-MAXDIST };
    double max[6] = {L/2.0, L/2.0, MAXDIST,MAXDIST,MAXDIST, MAXDIST };
    
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
   
    //// Imaginary
    iter=0;
    result_imag=0;
    helper.imaginary_part=true;
/*    do
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
    gsl_monte_miser_free(s);
    
    
    return 0;


}
