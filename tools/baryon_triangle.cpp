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
#include <vector>
#include <fstream>


gsl_rng* global_rng;
using namespace std;


double L = 5*Amplitude::FMGEV;
const double MAXDIST = 5 *Amplitude::FMGEV; // Max q distance from the center
const double MCINTACCURACY=0.1;
double probeB=4.0;
bool slopeb2=false;
bool elastic=false;
bool totxs = false;
bool normalize_area = true;
bool dipole = false; // compute dipole instead of baryon
const std::string datafile = "triangles";
struct inthelper
{
    IPGlasma *glasma;
    bool imaginary_part;
    std::vector<Vec> q1v; // probe quarks
    std::vector<Vec> q2v;
    std::vector<Vec> q3v;
};

//
// MC vector: (b1x,b1y,q1x,q1y,q2x,q2y)
// Location for vec q3 is fixed by q1+q2+q3=0
// q1,q2 are measured w.r.t. impact parameter
double inthelperf_mc(double* vec, size_t dim, void* p)
{
    inthelper *par = (inthelper*)p;
    
    Vec b(vec[0],vec[1]);
        // Generate random probe configuration
        // see https://math.stackexchange.com/questions/442418/random-generation-of-rotation-matrices

    /*    // First put quarks on an equilateral triangle
        Vec q1  (-probeB/2.0, 0., 0.);
        Vec q2  (probeB/2.0, 0., 0.);
        Vec q3  (0., sqrt(3.0/4.0)*probeB, 0.);

        // Shift center to origin
        Vec center = q1; center =center + q2; center = center+q3;
        center = center * (1.0/3.0);
        q1 = q1-center;
        q2=q2-center;
        q3=q3-center;

        if (dipole == true)
        {
            // Not a triangle but fixed distance probeB
            q1 = Vec(-probeB/2.0, 0, 0);
            q2 = Vec(probeB/2.0, 0, 0);
            // q3 not used
        }        

 // Sample random rotation vector
        double throt = std::acos(2.0*gsl_rng_uniform(global_rng)-1) ;
        double phi = gsl_rng_uniform(global_rng)*2.0*M_PI;

        Vec rotationvec (sin(throt)*cos(phi), sin(throt)*sin(phi), cos(throt));

        double rotate = gsl_rng_uniform(global_rng)*2.0*M_PI;
        q1.Rotate3D(rotationvec,rotate);
        q2.Rotate3D(rotationvec,rotate); 
        q3.Rotate3D(rotationvec,rotate);

        int probes=100;
*/
    double sum=0;
    int probes = par->q1v.size();
    probes=1000;
    for (int probe=0; probe<probes; probe++)
    {
        Vec q1 = par->q1v[probe]*probeB;
        Vec q2 = par->q2v[probe]*probeB;
        Vec q3 = par->q3v[probe]*probeB;
          
        double q1v[2]={q1.GetX()+b.GetX(), q1.GetY()+b.GetY()};
        double q2v[2]={q2.GetX()+b.GetX(), q2.GetY()+b.GetY()};
        double q3v[2]={q3.GetX() + b.GetX(), q3.GetY()+b.GetY()};
    
        complex<double> baryon;
        baryon = par->glasma->BaryonOperator(0.01, q1v,q2v,q3v);
        complex<double> res=baryon;
        if (totxs)
            res = 1.0 + baryon/6.0;
        if (elastic and totxs)
            res = (1.0 + baryon/6.0)*(1.0 + baryon/6.0);
        if (slopeb2 and totxs)
            res = b.LenSqr()/2.0 * (1.0 + baryon/6.0); 
        
    double normalization = std::pow(L, 2);
    if (normalize_area == false)    
        normalization=1;

    if (par->imaginary_part == false)
        sum += res.real()  / normalization;
    else
        sum += res.imag()  /normalization;
    }
    return sum/probes;
}



double inthelperf_dipole_mc(double* vec, size_t dim, void* p)
{
    inthelper *par = (inthelper*)p;
    
    Vec b(vec[0],vec[1]);
    Vec q1 = par->q1v[0];
    Vec q2 = par->q2v[0];
    cout << "I DON'T WORK!" << endl;
    cerr << "DIPOLE DOESNT WORK! " << endl;
    exit(1);
   
    double q1v[2]={q1.GetX()+b.GetX(), q1.GetY()+b.GetY()};
    double q2v[2]={q2.GetX()+b.GetX(), q2.GetY() + b.GetY()};
    
    double dipole = par->glasma->Amplitude(0.01, q1v, q2v);
    
    double normalization = std::pow(L, 2);
    if (normalize_area == false)    
        normalization=1;

    return dipole/normalization;

}
int main(int argc, char* argv[])
{
    int probes = 100;

    string fname = argv[1];
    int mcintpoints = StrToReal(argv[2]);
    L = StrToReal(argv[3])*Amplitude::FMGEV;
    string mode = argv[4];
    probeB=StrToReal(argv[5]);
    size_t dim = 2;
    if (mode == "baryon")
    {
        dipole = false;
        elastic=false;
        slopeb2=false;
    }
    else if (mode=="elastic_baryon")
    {
        dipole=false;
        elastic=true;
        slopeb2=false;
    }
    else if (mode == "elastic_slope_baryon")
    {
        dipole=false;
        elastic=true;
        slopeb2=true;
    }
    else if (mode == "dipole")
    {
        dipole = true;
        elastic=false;
    }
    else if (mode == "elastic_dipole")
    {
        dipole=true;
        elastic=true;
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


    // Read triangles
    std::vector<Vec> q1v;
    std::vector<Vec> q2v;
    std::vector<Vec> q3v;
    std::ifstream f(datafile.c_str());
    while (!f.eof() )
    {
        std::string line;
        std::getline(f, line);
        std::stringstream ss(line);
        if (line.length()<10) continue;
        double x,y,z;
        ss >> x; ss >> y; ss >> z;
        Vec q(x,y,z); q1v.push_back(q);

        ss >> x; ss >> y; ss >> z;
        Vec q2(x,y,z); q2v.push_back(q2);

        ss >> x; ss >> y; ss >> z;
        Vec q3(x,y,z); q3v.push_back(q3);

    }

    f.close();


    cout << "# Compute totxs: " << totxs << " normalize area: " << normalize_area << " periodic boundary conditions: " << periodicboundary << endl;
    cout << "# Baryon configs read: " << q1v.size() << endl;

    inthelper helper;
    helper.glasma = &glasma;
    helper.imaginary_part=false;
	

    helper.q1v=q1v; helper.q2v=q2v; helper.q3v=q3v; 

    gsl_monte_function fun;
    fun.params=&helper;
    double *min;
    double *max;
    if (dipole == false)
    {
        dim=2;
        min = new double[dim];
        max = new double[dim]; 
        min[0] = -L/2; min[1] = -L/2; // min[2] = -MAXDIST; min[3]=-MAXDIST; min[4]=-MAXDIST; min[5]=-MAXDIST;
        max[0] = L/2; max[1] = L/2; //max[2] = MAXDIST; max[3]=MAXDIST; max[4]=MAXDIST; max[5]=MAXDIST;

        fun.f = inthelperf_mc;
    }
    else
    {
        dim=2;
        min = new double[dim];
        max =new double[dim];
        min[0]=-L/2; min[1]=-L/2; // min[2] = -MAXDIST; min[3]=-MAXDIST;
        max[0] = L/2; max[1] = L/2; //max[2] = MAXDIST; max[3]=MAXDIST;
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
            }
         //gsl_monte_plain_integrate
             gsl_monte_miser_integrate
                (&fun, min, max, dim, mcintpoints, global_rng, s,
                               &result_real, &abserr);
            //cout << "relerr " << std::abs(abserr/result_real) << endl;
            
        } while (std::abs(abserr/result_real)>MCINTACCURACY);
    

    double real = result_real;
 
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
    cout <<real  << " " << 0 << endl;
    //gsl_monte_plain_free (s);
    delete [] min;
    delete[] max; 
    gsl_monte_miser_free(s);
    
    
    return 0;


}
