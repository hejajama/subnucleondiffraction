/* HERA sigmar fitter
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 *
 * Computes sigmar at given HERA kinematics
 *
 * 1st argument: JIMWLK solution dir
 * 2nd argument: step size
 * 3rd argument: max steps
 * 4th argument: HERA datafile
 *
 * HERA datafile with structure
 * Q^2 x inelasticity sigma_r error
 * Lines starting with # are ignored
 * inelasticity is the DIS y variable y=Q^2/(sx)
 *
 * Prints the computed reduced cross section values in the standard output, and the format is
 * Q^2 x y experimental_sigmar exp_error theory_sigmar
 */


#include "../src/dipole.hpp"
#include "../src/diffraction.hpp"


#include "../src/vector.hpp"
#include "../src/subnucleon_config.hpp"
#include "../src/ipglasma.hpp"
#include "../src/dis.hpp"
#include "../src/virtual_photon.hpp"

#include <gsl/gsl_rng.h>

#include <tools/tools.hpp>
#include <gsl/gsl_errno.h>

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sstream>
#include <gsl/gsl_errno.h>
using namespace std;
using namespace Amplitude;

double x0 = 0.01;
const double ds = 0.0004;
gsl_rng* global_rng;

double mean(vector<double> &v)
{
    double sum=0;
    for (unsigned int i=0; i<v.size(); i++) sum += v[i];
    return sum/v.size();
}

int main(int argc, char* argv[])
{
    MCINTPOINTS=1e5;
    gsl_set_error_handler(&ErrHandler);
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);
    cout << "# F2 fitter" << endl;
    if (argc < 6)
    {
        cout << "Syntax: " << argv[0] << " jimwlkdir step maxstep heradata alphas config" << endl;
        return 0;
    }
    
    string jimwlkdir = argv[1];
    int step = StrToInt(argv[2]);
    int maxstep = StrToInt(argv[3]);
    double alphas = StrToReal(argv[5]);
    string herafile = argv[4];
    double quarkmass = 1.4;
    int averages = StrToInt(argv[6]);
    
    bool fixed_coupling = false;
    double minx=0;
    if (alphas > 0)
    {
        // Fixed coupling
        fixed_coupling = true;
        double maxy =  M_PI*M_PI / alphas * maxstep * ds;
        minx = x0*exp(-maxy);
        cout << "# x range in fit: " << x0 << " - " << minx << ", # of configs" << averages << endl;
    }
    
    bool scale_x = false;   // scale bjorken x to take into account the quark mass
    bool include_light = true;
    
    int points=0;
    double chisqr=0;
    
    
    if (scale_x)
        cout <<"# Shifting Bjorken-x to x*(1 + 4mf^2/Q^2)" << endl;
    
    // Rear HERA data
    vector<double> xvals, yvals, qsqrvals, expvals, experrors;
    //NB: yvals is now vector of the inelasticities, y=Q^2/(sx)
    ifstream file(herafile.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << file << endl;
        return -1;
    }
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#")
        continue;
        string x,qsqr,y,sigmar,err;
        stringstream l(line);
        l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
        //if (std::abs(StrToReal(f2)) <1e-10) continue;
        if (StrToReal(x)>x0 or StrToReal(x)<minx) continue;
        qsqrvals.push_back(StrToReal(qsqr)); xvals.push_back(StrToReal(x));
        yvals.push_back(StrToReal(y));
        expvals.push_back(StrToReal(sigmar)); experrors.push_back(StrToReal(err));
    }
    file.close();
    
    
    
    
    cout << "# Q^2 [GeV^2]  x  y  HERA-\\sigma_r  HERA-err theory-\\sigma_r(light c b) " << endl;
    // Compute reduced cross section
    for (unsigned int i=0; i<xvals.size(); i++)
    {
        double x = xvals[i];
        double sqrts = std::sqrt( qsqrvals[i]/(x * yvals[i]) );
        
        if (scale_x)
            x = x   * (1.0 + 4.0*SQR(quarkmass)/qsqrvals[i]);
        
        // Calculate reduced xs below and above, and interpolate
        // Fixed coupling: step = alphas*ln(x0/x) / (pi^2 * ds)
        int steps_upper = 0;
        int steps_lower = 0;
        int evolsteps;
        
        if (fixed_coupling)
        {
            evolsteps = alphas * log(x0/x)/ (M_PI*M_PI*ds);
            steps_lower = evolsteps - evolsteps%step;
            steps_upper = evolsteps + (step - evolsteps%step);
            //cout << " Evolsteps " << evolsteps << " lower " << steps_lower << " upper " << steps_upper << endl;
        }
        
        VirtualPhoton photon;
        photon.SetQuark(Amplitude::C, quarkmass);
        if (include_light)
            photon.AddQuark(Amplitude::LIGHT, 0.14);
        vector<double> xs_t_upper;
        vector<double> xs_l_upper;
        vector<double> xs_t_lower;
        vector<double> xs_l_lower;
        //for(int conf = 0; conf < averages; conf++)
        int conf = averages; // Do only single configuration, average later separately!
        {
            stringstream fname_upper;
            fname_upper << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_upper;
            stringstream fname_lower;
            fname_lower << jimwlkdir << "/V-" << conf << ".txt_steps_" << steps_lower;
            //cout << fname_upper.str() << " " << fname_lower.str() << endl;
            
            IPGlasma dipole_upper(fname_upper.str());
            IPGlasma dipole_lower(fname_lower.str());
            
            Diffraction f2upper(dipole_upper, photon);
            Diffraction f2lower(dipole_lower, photon);
            
            double t_up = 4.0*M_PI*f2upper.ScatteringAmplitude(x, qsqrvals[i], 0, T);
            double l_up = 4.0*M_PI*f2upper.ScatteringAmplitude(x, qsqrvals[i], 0, L);
            xs_t_upper.push_back(t_up);
            xs_l_upper.push_back(l_up);
            
            double t_low = 4.0*M_PI*f2lower.ScatteringAmplitude(x, qsqrvals[i], 0, T);
            double l_low = 4.0*M_PI*f2lower.ScatteringAmplitude(x, qsqrvals[i], 0, L);
            xs_t_lower.push_back(t_low);
            xs_l_lower.push_back(l_low);
            
          

            
        }
        
        double xs_t_u =mean(xs_t_upper);
        double xs_l_u = mean(xs_l_upper);
        double f2_upper = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_u + xs_l_u);
        double fl_upper = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_u );
        double sigmar_upper = f2_upper - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_upper;
        
        double xs_t_l =mean(xs_t_lower);
        double xs_l_l = mean(xs_l_lower);
        double f2_lower = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * (xs_t_l + xs_l_l);
        double fl_lower = qsqrvals[i]/(4.0*M_PI*M_PI*ALPHA_e) * ( xs_l_l );
        double sigmar_lower = f2_lower - yvals[i]*yvals[i]/(1.0 + pow(1.0-yvals[i],2.0))*fl_upper;
        
        // Interpolate in s ~ ln 1/x
        double sigmar = sigmar_lower + (evolsteps-steps_lower) / (steps_upper - steps_lower) * (sigmar_upper - sigmar_lower);
        
        cout << qsqrvals[i] << " " << xvals[i] << " " << yvals[i] << " " << expvals[i] << " " << experrors[i] << " " << sigmar << endl;
        
        points++;
        chisqr += pow((sigmar - expvals[i])/experrors[i],2.0);
    }
    
    cout << "#Chi^2/N " << chisqr/points << " points " << points << endl;
    
}
