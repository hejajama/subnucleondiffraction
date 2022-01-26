/*
 * Nucleus that consists of nucleons
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "nucleons.hpp"
#include "vector.hpp"
#include "ipsat_proton.hpp"
#include "subnucleon_config.hpp"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <sstream>
#include <string>
#include "interpolation.hpp"
#include <fstream>

#define LINEINFO __FILE__ << ":" << __LINE__

using namespace std;




double Nucleons::Amplitude(double xpom, double q1[2], double q2[2] )
{
    // Calculate scattering amplitude in Glauber approach
    // Idea: the S matrix is a product of dipole-nucleon S matrixes
    // S(r,b) = prod_i S_p(r, b-b_i);
    // where b_i is the center of nucleon i
    
    double smat=1;
    Vec qv1(q1[0],q1[1]);
    Vec qv2(q2[0],q2[1]);
    Vec b = qv1 + qv2;
    b = b * 0.5;
    Vec r= qv2 - qv1;
    if (A==2 and deuteron_structure == TUBE)
    {
        double tp_deuteron = DeuteronTubeDensity(b);
        return ((Ipsat_Proton*)nucleons[0])->Amplitude_Tp(xpom, r.Len(), tp_deuteron );
    }
    for (int i=0; i<A; i++)
    {
        // Calculate the quark and antiquark coordinates in the frame where the nucleon i is at the origin
        Vec new_q1  = qv1 - nucleon_positions[i];
        Vec new_q2 = qv2 - nucleon_positions[i];
        smat = smat * (1.0 - nucleons[i]->Amplitude(xpom, new_q1, new_q2));
        
    }
    return 1.0-smat;
}

// Helper to read VMC density function
Interpolator ReadVMC(std::string fname="h2.density");

void Nucleons::InitializeTarget()
{
    nucleon_positions.clear();
    
    // Handle special cases
    if (A==2)
    {
        double maxr = 25*FMGEV;
        
        nucleons[0]->InitializeTarget();
        nucleons[1]->InitializeTarget();
        // Sample difference
        
        // Woods Saxon separately
        if (DeuteronWF == WoodsSaxon)
        {
            Vec tmp;
            do {
                Vec tmpvec (2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr);
                tmp=tmpvec;
            } while (gsl_rng_uniform(global_rng) > WS_unnorm(tmp.Len())); // WS distribution!
            nucleon_positions.push_back(tmp);
            Vec tmp2 = tmp*(-1.0);
            nucleon_positions.push_back(tmp2);
            
        }
        
        else    // Hulthen, or VMC, implemented in DeuteronWaveFunction
        {
            double probability=0;
            Vec tmp;
            do {
                
                Vec tmpvec (2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                            2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr);
                tmp=tmpvec;
                probability = DeuteronWaveFunction(tmp.Len());
                
            } while (gsl_rng_uniform(global_rng) > probability);
            // Now vec is from proton to neutron, put origin at the center
            Vec proton = tmp*0.5;
            Vec neutron = proton*(-1);
            //proton.SetZ(0);
            //neutron.SetZ(0);
            nucleon_positions.push_back(proton);
            nucleon_positions.push_back(neutron);
        
        }

    }
    // He 3
    else if (A==3)
    {
        nucleons[0]->InitializeTarget();
        nucleons[1]->InitializeTarget();
        nucleons[2]->InitializeTarget();
        // Nucleon coordinates are in file he3.dat
        // Format: x y z x y z x y z
        std::ifstream f("he3.dat");
        
        if (!f.is_open())
        {
            std::cerr << "Could not open file he3.dat " << std::endl;;
            exit(1);
            return;
        }
        std::string line;
        
        int index=-1;
        bool found=false;
        while(!f.eof() )
        {
            index++;
            std::getline(f, line);
            if (index != he3_id)
                continue;
            
            // This is the correct line
            found=true;
            double x1,y1,z1,x2,y2,z2,x3,y3,z3;
            stringstream ss(line);
            ss>>x1; ss>>y1; ss>>z1;
            ss>>x2; ss>>y2; ss>>z2;
            ss>>x3; ss>>y3; ss>>z3;
            Vec n1(x1*FMGEV,y1*FMGEV); Vec n2(x2*FMGEV,y2*FMGEV); Vec n3(x3*FMGEV,y3*FMGEV);
            nucleon_positions.push_back(n1);
            nucleon_positions.push_back(n2);
            nucleon_positions.push_back(n3);
            break;
            
        }
        if (found==false)
        {
            cerr << "Did not found He3 config " << he3_id << endl;
            exit(1);
        }
    }
    else
    {
        for (int i=0; i<A; i++)
        {   
            nucleons[i]->InitializeTarget();
            double maxr = 3.0*ws_ra;
            
            Vec tmp;
            do {
                Vec tmpvec (2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                                2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr,
                                2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr);
                tmp=tmpvec;
            } while (gsl_rng_uniform(global_rng) > WS_unnorm(tmp.Len())); // WS distribution!
            nucleon_positions.push_back(tmp);
        }
    }
    
}

Nucleons::Nucleons(std::vector<DipoleAmplitude*> nucleons_) : VMC_interpolator(ReadVMC("h2.density"))
{
    A=nucleons_.size();
    cout << nucleons_[0]->InfoStr();
    nucleons=nucleons_;
    //DeuteronWF = Hutlhen;
    DeuteronWF = VMC;
    
    if (A > 2)
    {
        ws_delta=0.54*FMGEV;
        ws_ra = 1.12 * std::pow(A, 1.0/3.0) * FMGEV;
    }
    else
    {
        // PHOBOS parameters
        ws_ra = 0.01 * FMGEV;
        ws_delta = 0.5882 * FMGEV;
    }
    he3_id=-1;
    deuteron_structure = NUCLEONS;

}


double Nucleons::WS_unnorm(double r )
{
    return 1.0 / (1+exp((r-ws_ra)/ws_delta));
    
}

std::string Nucleons::InfoStr()
{
    std::stringstream ss;
    
    if (A==2)
    {
        ss << "#DipoleAmplitude: Deuteron nucleus, nucleon coordinates (" << nucleon_positions[0].GetX() << ", " << nucleon_positions[0].GetY() << ") and (" <<  nucleon_positions[1].GetX() << ", " << nucleon_positions[1].GetY() << endl;
        ss << nucleons[0]->InfoStr() << endl << nucleons[1]->InfoStr() << endl;
        if (DeuteronWF == Hulthen)
            ss << "# Deuteron wave function: Hulthen" << endl;
        else if (DeuteronWF == ExtendedHulthen)
            ss << "# Deuteron wave function: ExtendedHulthen" << endl;
        else if (DeuteronWF == WoodsSaxon)
            ss << "# Deuteron wave function:  WoodsSaxon" << endl;
        else if (DeuteronWF == VMC)
            ss << "# Deuteron wave function: VMC" << endl;
        
        if (deuteron_structure == NUCLEONS)
            ss << "# Deuteron = proton + neutron " << endl;
        else
            ss << "# Deuteron = tube connecting p and n" << endl;
    }
    else if (A==3)
    {
        ss << "#DipoleAmplitude: He3 configuration " << he3_id << endl;
        ss << "# (" << nucleon_positions[0].GetX() << ", " << nucleon_positions[0].GetY() << "), (" <<nucleon_positions[1].GetX() << ", " << nucleon_positions[1].GetY() << "), (" << nucleon_positions[2].GetX() << ", " << nucleon_positions[2].GetY() << ")" << endl;
     }
    else
    {
        ss << "#DipoleAmplitude: Nucleus cosisting of " << A << " nucleons, nucleon 0 info: " << endl;
        ss << nucleons[0]->InfoStr() << endl;
    }
    return ss.str();
}

Nucleons::~Nucleons()
{
    // Nucleons are allocated in main.cpp but should be freed here
    // If DGLAPdist was allocated in main.cpp, it is also freed there
    for (int i=0; i<A; i++)
        delete nucleons[i];
}

/*
 * Deuteron wave function, r is a 3d vector
 * used in IP-Glasma calculations in 1304.3403, 
 * orig. ref. Ann. Rev. Nucl. Part. Sci. 57, 205 (2007).
 * Probability is wavef^2
 *
 * Also supports Variational Monte Carlo (VMC) from https://www.phy.anl.gov/theory/research/density2/
 */


double Nucleons::DeuteronWaveFunction(double r)
{
    if (r<1e-10)
        return 0;   // Forbid zero distance
    
    if (DeuteronWF == Hulthen)
    {
        double a = 0.228/FMGEV; // 0.228 1/fm = 0.228/5.068 GeV
        double b =1.18/FMGEV;
        double probability_amp = 1.0/sqrt(2.0*M_PI) * sqrt(a*b*(a+b))/(b-a) * (exp(-a*r) - exp(-b*r))/r;
        
        if (1000.0*probability_amp*probability_amp > 1)
        {
            cerr << "DeuteronWF rejection sampling gets prob.amp. " <<1000.0*probability_amp*probability_amp << endl;
            exit(1);
        }
        return 1000.0*probability_amp*probability_amp;
        // Todo: maximum of this is very small, so rejection sampling does not work well without prefactor
    }
    else if (DeuteronWF == ExtendedHulthen)
    {
        // This parametrization is already squared, and assumes r is in fm
        r = r / FMGEV;
        if (r<0.02)
            return 0;
        //Phys. Rev. 151, 772
        double alpha = 0.2338;
        double Cvals[4] = {-0.63608, -6.6150, 15.2162, -8.9651};
        double evals[4] = {5.733*alpha, 12.844*alpha, 17.331*alpha, 19.643*alpha};
        
        double sum=0;
        for (int i=0; i<4; i++)
            sum += Cvals[i]*exp(-evals[i]*r);
        sum += exp(-alpha*r);
        
        return sum*sum /(r*r);   // This is always <1, so rejection sampling works
        // Note that here the parametrization in the Ref. is for the r distribution and thus includes Jacobian, so we have to remove it by normalizing by r^2
        // I have checked that at r>0.02fm this function is <1, thus rejection sampling works
        // sum*sum, as we return probability which is wavef^2
    }
    else if (DeuteronWF == VMC)
    {
        
        if (r > VMC_interpolator.MaxX())
        {
            return 0;
        }
        
        // NOTE Vmc data is scaled such that interpolator takes p-n-distance
        double rhorp = VMC_interpolator.Evaluate(r);

        if (rhorp < 0)
        {
            if (r/FMGEV < 20)
                cout << "Skip r that gives interp=" << rhorp << ", this is r = " << r/FMGEV << " fm " << endl;
            return 0;    // Crazy things happen at very large r, in practice at r >14fm
        }
        
        
        if (rhorp > 1)
        {
            cerr << "Crazy rho_rp(r=" << r << ") = " << rhorp << " at " << LINEINFO << endl;
            exit(1);
        }

        return rhorp;   // As far as I understand, there is no need for Jacobian here
        
    }
    
    std::cerr << "Unknown deuteron wave function!" << LINEINFO << endl;
    exit(1);
}

struct inthelper_deuteron_z
{
    Vec b;
    Vec nucleon1;
    Vec nucleon2;
    double B_p;
};
double inthelperf_deuteron_z(double z, void* p)
{
    inthelper_deuteron_z* par = (inthelper_deuteron_z*)p;
    par->b.SetZ(z);
    
    // Distance from the line connecting quarks
    Vec q1_to_q2 = par->nucleon2- par->nucleon1;
    Vec q2_to_q1 = par->nucleon1 - par->nucleon2;
    Vec q1_to_point = par->b - par->nucleon1;
    Vec q2_to_point = par->b - par->nucleon2;
    // Test if we are in between
    // Study projections q1->point on direction q1->q2, must be positive
    // similarly q2->point projected on q2->q1 must be positive
    double proj1 = q1_to_point*q1_to_q2;
    double proj2 = q2_to_point*q2_to_q1;
    if (proj1>0 and proj2>0)
    {
        double dsqr_line = ( q1_to_point.LenSqr() * q1_to_q2.LenSqr() -pow( q1_to_point*q1_to_q2, 2.0 )) / q1_to_q2.LenSqr() ;
        double density = std::exp(-dsqr_line / (2.0*par->B_p)) * 1.0 / (2.0*std::sqrt(2.0) * pow(par->B_p * M_PI, 3.0/2.0));
        return density;
    }
    else
    {
        double dsqr = std::min( q1_to_point.LenSqr(), q2_to_point.LenSqr());
        return std::exp(-dsqr / (2.0*par->B_p)) * 1.0 / (2.0*std::sqrt(2.0) * pow(par->B_p * M_PI, 3.0/2.0));
    }
    
    
}
double Nucleons::DeuteronTubeDensity(Vec b)
{
    // Connect p and n by a Gaussian tube
    
    // small optimization if we are faaar away
    Vec d1 = b - nucleon_positions[0]; d1.SetZ(0);
    Vec d2 = b - nucleon_positions[1]; d2.SetZ(0);
    if (d1.LenSqr() > 30*30 and d2.LenSqr() > 30*30)
        return 0;
    
    
        inthelper_deuteron_z helper;
        helper.b = b;
        helper.B_p=4;
        helper.nucleon1  = Vec( nucleon_positions[0].GetX(), nucleon_positions[0].GetY(), nucleon_positions[0].GetZ());
        helper.nucleon2  = Vec( nucleon_positions[1].GetX(), nucleon_positions[1].GetY(), nucleon_positions[1].GetZ());
        const int INTPOINTS_ZINT = 5;
        gsl_function f;
        f.function = &inthelperf_deuteron_z;
        f.params = &helper;
        double tuberesult,error;
        gsl_integration_workspace * w =  gsl_integration_workspace_alloc (INTPOINTS_ZINT);
        int status = gsl_integration_qag (&f, -40, 40, 0, 0.01, INTPOINTS_ZINT, GSL_INTEG_GAUSS15,
                             w, &tuberesult, &error);
        gsl_integration_workspace_free(w);
        // return with some reasonable normalization
    
        if (status)
            cerr << "Tube z integral failed, res " << tuberesult << " relerr " << error/tuberesult << endl;
    
	// best normalization is 0.46???
        return tuberesult*0.43;
    
    //return std::max(tuberesult, endresult);
    
}

void Nucleons::SetHeId(int i)
{
    if (i<0 or i>13698)
    {
        cerr << "He3 id " << i << " can not be used! Must be between 0 and 13698" << endl;
        exit(1);
    }
    he3_id=i;
}

std::vector<DipoleAmplitude*> Nucleons::GetNucleons()
{
    return nucleons;
}


Interpolator ReadVMC(std::string fname)
{
    std::ifstream f(fname.c_str());
    if (!f.is_open())
    {
        std::cerr << "Could not open file " << fname << " " << LINEINFO << std::endl;;
        exit(-1);
    }
    
    std::string line;
    std::vector<double> rvals;
    std::vector<double> rhovals;
    
    // Interpolate to 0
    rvals.push_back(0);
    rhovals.push_back(0);
    
    while(!f.eof() )
    {
        std::getline(f, line);
        if (line[0]=='#' or line.length() < 10)   // Comment line, or empty
        {
            continue;
        }
        
        // Parse
        std::stringstream ss(line);
        // File syntax
        // R    INTERP(R,RHOL,RLH)    RHORP      DRHORP
        double r,tmp,rho;
        ss >> r;
        ss >> tmp;
        ss >> rho;
        rvals.push_back(2.0*r*FMGEV); // Scale, so this is p-n-distance
        rhovals.push_back(rho);
        
    }
    
    f.close();
    
    Interpolator interp(rvals,rhovals);
    return interp;
}

