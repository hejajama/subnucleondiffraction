/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a proton that consists of quarks
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#include "ipsat_proton.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <tools/tools.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include "subnucleon_config.hpp"
#include "mz_ipsat/dipoleamplitude.hpp"

// IPsat 2012
#ifdef USE_FORTRAN_IPSAT12
extern "C" {
       double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
     };
#endif

int IPSAT12_PAR = 2;    // m_c=1.4 GeV

using std::cout; using std::endl;

const double FMGEV = 5.06778;
double MAXR_SKEW = 20;  // Dont calculate skew at larger r, as it didnt work

const int INTPOINTS_ZINT = 7;   // Depth of z integral subintervals when projecting exponential
// distribution into 2d


double Ipsat_Proton::Amplitude( double xpom, double q1[2], double q2[2])
{
    // quark transverse coodinates are now nucleons[i].GetX() and GetY()
    Vec q(q1[0], q1[1]);
    Vec qbar(q2[0],q2[1]);
    
    Vec b = q + qbar;
    b = b*0.5;
    
    Vec r = q - qbar;
    
    double tpsum = 0;
    tpsum = Density(b);
    double blen = b.Len();
    // pi^2/2Nc as * xg
    double ipsat_exponent = 0;
    if (ipsat == IPSAT06)
    {
        ipsat_exponent = gdist->Gluedist(xpom, r.LenSqr());
    }
    else if (ipsat == IPSAT12)
    {
#ifndef USE_FORTRAN_IPSAT12
        cerr << "IPsat12 support is not complied (requires Fortran compiler!)" << endl;
        exit(1);
#else
        // dipole_amplitude(xBj, r, b, parameterSet) gives amplitude 2(1 - exp(c*T(p)))
        double tmpb=0;  double tmpr = r.Len();
       	// If we do not have fluctuations, do not modify geometry
        if (B_q==4.0 and B_p == 0 and saturation  and Qs_fluctuation_sigma==0)
            return dipole_amplitude_(&xpom, &tmpr, &blen, &IPSAT12_PAR)/2.0;
 
        // We have to calculate "gluedist" as in case of ipsat06
        // par 1: m_c=1.27,   2: m_c=1.4
        double n = dipole_amplitude_(&xpom, &tmpr, &tmpb, &IPSAT12_PAR)/2.0;
       
 
        double c = -std::log(1.0-n);
        
        if (std::isnan(c) or std::isinf(c))
        {
            // We have so large r, that basically n=1 and c blows up, these should not matter
            // as wave function cuts these out anyway, but we can just set amplitude to 1
            //
            // NOTE: When you calculate F_2, be very careful! This may have significant effect
            // Also, in case of light mesons, there could be some effect!
            
            // Fall back to round proton! This gives approximately 1 if impact parameter is not crazy
            double tmpb = b.Len();
            return dipole_amplitude_(&xpom, &tmpr, &tmpb, &IPSAT12_PAR)/2.0;
            //return 1.0;
        }
        
        double tp = 1.0/(2.0*M_PI*4.0)*std::exp(- tmpb*tmpb / (2.0*4.0));
        c /= tp;
        
        ipsat_exponent = c / r.LenSqr();
#endif
    }
    else if (ipsat == MZSAT or ipsat == MZNONSAT)
    {
        const double Nc=3;
        ipsat_exponent = M_PI * M_PI / (2.0*Nc) * mzipsat->Alphas_xg(xpom, mzipsat->MuSqr(r.Len()));
    }
#ifdef USE_LCPT_DIPOLE
    else if (ipsat == LCPT)
    {
		double phirb = std::acos(r*b/r.Len()*b.Len());
        double n = lcpt_dipole->Evaluate(r.Len(), b.Len(),phirb);
        return n; // Note: does not support geometry params
    }
#endif
    else
    {
        std::cerr << "UNKNOWN IPSAT VERSION!" << std::endl;
        exit(1);
    }
    
    double skew = 1.0;
    if (skewedness and r.Len() < MAXR_SKEW)
    {
        double skew_lambda = LogDerivative_xg(xpom, r.Len());
        if (isnan(skew_lambda))
            cerr << "NaN skew, probably too large r=" << r.Len() << "? Neglecting skew now.. " << LINEINFO << endl;
        else
            skew = Skewedness(skew_lambda);
    }
    double tmpb=b.Len(); double tmpr=r.Len();
    
    
    if (saturation)
        return 1.0 - std::exp( - r.LenSqr() * 1.0/quarks.size() * GetQsFluctuation(b.GetX(),b.GetY()) * skew* ipsat_exponent * tpsum  );
    else
        return r.LenSqr() * 1.0/quarks.size() * GetQsFluctuation(b.GetX(),b.GetY()) * skew* ipsat_exponent * tpsum;
    

   
    return 0;
    
}

/*
 * Proton density
 * T_p(b), or something that replaces it
 */
double Ipsat_Proton::Density(Vec b)
{
    double density=0;
    if (proton_structure == QUARKS)
    {
        // Need to calculate sum_i T_p(|b-b_i|), where b is the center of the dipole
        // and b_i is the center of the quark i
        for (unsigned int i=0; i < quarks.size(); i++)
        {
            Vec projection (quarks[i].GetX(), quarks[i].GetY());
            Vec deltab = b - projection;
            density = density + QuarkThickness(deltab.Len(), i);
        }
        return density;
    } else if (proton_structure == CENTER_TUBES)
    {
        density = FluxTubeThickness(b);
        density *= quarks.size(); // Because this sum is later normalized by 1/N_q
        return density;
        
    }
    else{
        cerr << "Unkown proton profile!" << endl;
        exit(1);
    }
    
}

double Ipsat_Proton::Amplitude_Tp(double xpom, double r, double tp)  // Dipole ampliutde at fixed tp
{
    if (ipsat == IPSAT12)
    {
#ifndef USE_FORTRAN_IPSAT12
        cerr << "IPsat12 support is not complied (requires Fortran compiler!)" << endl;
        exit(1);
#else
        if (!saturation)
        {
            std::cerr << "Nonsat is not defined for ipsat2012" << std::endl;
            exit(1);
        }
        // dipole_amplitude(xBj, r, b, parameterSet) gives amplitude 2(1 - exp(c*T(p)))
        // We have to calculate "gluedist" as in case of ipsat06
        double tmpb=0;
        // par 1: m_c=1.27,   2: m_c=1.4
        double n = dipole_amplitude_(&xpom, &r, &tmpb, &IPSAT12_PAR)/2.0;
        
        double c = std::log(1.0-n);
        
        if (std::isnan(c) or std::isinf(c))
        {
            // We have so large r, that basically n=1 and c blows up, these should not matter
            // as wave function cuts these out anyway, but we can just set amplitude to 1
            return 1.0;
        }
        
        double tp0 = 1.0/(2.0*M_PI*4.0)*std::exp(- tmpb*tmpb / (2.0*4.0));
        c /= tp0;
        c *= tp;
        
        return 1.0 - std::exp(c);
#endif
    }
    else
    {
        cerr << "Ipsat_Proton::Amplitude_Tp is only defined for IPSAT12" << endl;
        exit(1);
    }
    
    return -1;
    
}


void Ipsat_Proton::InitializeTarget()
{
    double smallestdist=999; double largestdist=0;
    double tmpdist;

    
    
    quarks.clear();
    quark_bp.clear();
    
    // Correlated sampling
    if (shape == ALBACETE)
    {
        /*
         * Implement arXiv:1605.09176.
         * Quark positions s_i are sampled from a distribution
         * prod_i<j, i,j=1...n_q  (1-exp(-mu*|s_i - s_j|^2/R^2))
         * Where s_i are sampled also from a Gaussian distribution that has width R
         * Note that in this case we keep the proton center-of-mass at the origin, thus
         * we only sample two quarks which then fixes the position of the third quark
         * Implemented using simple rejection method
         *
         * Note that we use a Gaussian distribution exp(-b^2/2B) when sampling the
         * constituent quark radii, which is different parametrization form Ref.
         */
        double r_c = 0.3 * 5.068;   // Repulsive core radius 0.3fm
        if (B_p < 1e-5)
        {
            std::cerr << "Does not support zero size proton!" << std::endl;
            exit(1);
        }
        else
        {
            std::vector<Vec> tmpquarkvec;
            int rejections=-1;
            do{
                rejections++;
                tmpquarkvec.clear();
                double xsum=0;
                double ysum=0;
                for (int i=0; i<number_of_quarks-1; i++)
                {
                    // Sample quark i
                    double x,y;
                    x=gsl_ran_gaussian(global_rng, std::sqrt(B_p));
                    y=gsl_ran_gaussian(global_rng, std::sqrt(B_p));
                    xsum +=x;
                    ysum +=y;
                    Vec s(x,y);
                    tmpquarkvec.push_back(s);
                }
                // Calculate last quark position by fixing origin=center of mass
                Vec q( -xsum, -ysum );
                tmpquarkvec.push_back(q);
                
            } while (gsl_rng_uniform(global_rng) > RepulsiveGaussianDistribution(tmpquarkvec, r_c));
            cout << "# Rejected " << rejections << " configurations due to short-range repulsive correlations" << endl;
            // Save accepted quarks
             for (int i=0; i<tmpquarkvec.size(); i++)
            {
                quark_bp.push_back(B_q);
                quarks.push_back(tmpquarkvec[i]);
                quarks3d.push_back(tmpquarkvec[i]);
            }
            
        }
     }
    else if (shape == MORELAND)
    {
        // Scott's Python script converted to C++
        // see original in scott_sampler/sampler.py
        // Numbers are in fm, converted to GeV^-1 at the end
        // NOTE: Some parameters also copied to SampleQSFluctuations
        double radius = 0.87;
        double structure_parameter = 0.3;
        unsigned int constituent_number=5;
        double sigma_fluct=0.81;
        double constituent_width = 0.2 + structure_parameter*(radius - 0.2);
        double sampling_radius = std::sqrt( std::pow(radius,2) - std::pow(constituent_width,2));
        double gamma_shape = 1./(constituent_number*std::pow(sigma_fluct,2) );
        
        for (unsigned int i=0; i < constituent_number; i++)
        {
            double qx = gsl_ran_gaussian(global_rng, sampling_radius);
            double qy = gsl_ran_gaussian(global_rng, sampling_radius);
            Vec tmpvec(qx*FMGEV,qy*FMGEV,0);
            Vec tmpvec3d(qx*FMGEV,qy*FMGEV,0);  // Not really 3d!
            quarks.push_back(tmpvec);
            quark_bp.push_back(constituent_width*constituent_width*FMGEV*FMGEV);
            quarks3d.push_back(tmpvec3d);
            
        }
    }
        // Sample uncorrelated quark positions
     else {
         
        // Sample 3 quarks
        for (int i=0; i<3; i++)
        {
            // Radius from uniform distribution
            double radius=0;
            double maxr = 30;
            
            if (shape == GAUSSIAN)
            {
                double x,y,z;
                if (B_p < 1e-5)
                {
                    x=y=z=0;
                    //radius=0;
                }
                
                else
                {
                    x=gsl_ran_gaussian(global_rng, std::sqrt(B_p));
                    y=gsl_ran_gaussian(global_rng, std::sqrt(B_p));
                    z=gsl_ran_gaussian(global_rng, std::sqrt(B_p));
                    /*
                    do{
                        radius = gsl_rng_uniform(global_rng) * maxr;
                    } while (gsl_rng_uniform(global_rng) > GaussianRadiusDistribution(radius));
                     */
                }
                // Sample angle
                //double angle = 2.0*M_PI*gsl_rng_uniform(global_rng);
                //Vec tmpvec(radius*std::cos(angle), radius*std::sin(angle));
                Vec tmpvec(x,y,0);
                Vec tmpvec3d(x,y,z);
                quarks.push_back(tmpvec);
                quarks3d.push_back(tmpvec3d);
                quark_bp.push_back(B_q);
            }
            else if (shape == EXPONENTIAL)
            {
                // We have to sample x,y,z separately
                double x,y,z;
                if (B_p < 1e-5)
                {
                    x=y=z=0;
                }
                else
                {
                    do{
                        x = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
                        y = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
                        z = 2.0*(gsl_rng_uniform(global_rng)-0.5)*maxr;
                    } while (gsl_rng_uniform(global_rng) > ExponentialDistribution(x,y,z));
                }
                Vec tmpvec(x,y,0);
                quarks.push_back(tmpvec);
                Vec tmpvec3d(x,y,z);
                quarks3d.push_back(tmpvec3d);
                quark_bp.push_back(B_q);
            }
            
        }
     }
    
    SampleQsFluctuations();

    // set center of mass to origin
    if (origin_at_center_of_mass)
    {
        Vec centerofmass;
        Vec centerofmass3d;
        for (int i=0; i<quarks.size(); i++)
        {
            centerofmass += quarks[i];
            centerofmass3d+=quarks3d[i];
        }
        centerofmass*=1.0/quarks.size();
        centerofmass3d*=1.0/quarks.size();
        for (int i=0; i<quarks.size(); i++)
        {
            quarks[i] -= centerofmass;
            quarks3d[i] -= centerofmass3d;
        }
    }
    
    
    // Calculate center of the quark triangle
    if (proton_structure == CENTER_TUBES)
    {
        center = GeometricMedian(quarks);
        center3d = GeometricMedian(quarks3d);
        NormalizeFluxTubeThickness();
    }
    
}

/*
 * Sample Q_s fluctuations at each point in transverse plane
 */
void Ipsat_Proton::SampleQsFluctuations()
{
    qs_fluctuation_coordinates.clear();
    qs_fluctuation.clear();
    qs_fluctuations_quarks.clear();
    
    if (std::abs(Qs_fluctuation_sigma) < 1e-10)
    {
        for (int i=0; i < quarks.size(); i++)
            qs_fluctuations_quarks.push_back(1.0);
        return; // No fluctuations
    }
    
    // ln Q_s^2 fluctuates according to a log-normal distribution
    // 1/(sqrt(2pi) sigma) exp[ -ln^2 x / (2 sigma^2) ],
    // where x = Q_s / <Q_s>
    
    // As we want that the average Q_s does not change, the sampled factor
    // is normalized by the arithmetic mean of the log-normal
    // distribution
    // The arithmetic mean of log0noral distribution exp(- ln^2 x/(2\sigma^2))
    // is exp(1/2 * sigma^2), and in our case sigma=0.5
    // This gives 1.133148
    
    double lognormal_mean = std::exp(0.5*Qs_fluctuation_sigma*Qs_fluctuation_sigma);;
    
    if (fluctuation_shape == FLUCTUATE_QUARKS)
    {
        
        
    
        // Saturation scale of each quark fluctuates from a Gaussian distribution
        // Note that as Q_s^2 ~ xg * T_q, and we get ln Q_s^2 fluctuations from
        // a Gaussian distribution, this same distribution can be used to describe
        // the normalization fluctuations of quarks
        
        
        int nq  = quarks.size();
        double sum=0;
        for (int i=0; i<nq; i++)
        {
            double f = 1;
            if (shape == MORELAND)
            {
                double sigma_fluct=0.81;
                double gamma_shape = 1./(nq*std::pow(sigma_fluct,2) );
                f=gsl_ran_gamma(global_rng, gamma_shape, 1./gamma_shape);
            }
            else
            {
                f = gsl_ran_gaussian(global_rng, Qs_fluctuation_sigma);
                f = std::exp(f)/lognormal_mean;
            }
            qs_fluctuations_quarks.push_back(f);
            sum+=f;
        }
        cout << "# Sampled " << nq << " quark fluctuations ";
        for (int i=0; i<nq; i++) cout << qs_fluctuations_quarks[i] << " ";
        cout << " average Q_s^2 fluctuation " << sum/nq << endl;
        
    }
    else if (fluctuation_shape == LOCAL_FLUCTUATIONS)
    {
        // Note for small grids: The saturation scale at point (x,y) is found by
        // seeking from the coordinates lists the index for which coordinates[i]<x<coordinates[i+1]
        // Thus, in practice having points=3 corresponds to 2x2 grid, as the leftmost values are never used
        double size = 10;
        int points = 11;    // x*x grid
        double step = (2.0*size)/(points-1);
        for (double x=-size; x<=size; x+=step)
            qs_fluctuation_coordinates.push_back(x);
        
        
        //double fluct=gsl_ran_gaussian(global_rng, Qs_fluctuation_sigma);
        double fluct_coef_sum=0; int pts=0; double stdev = 0;
        for (int yind=0; yind<points; yind++)
        {
            std::vector<double> row;
            for (int xind=0; xind<points; xind++)
            {
                double fluct;
                if (Qs_fluctuation_sigma > 0)
                {
                    double f = gsl_ran_gaussian(global_rng, Qs_fluctuation_sigma);
                    fluct = std::exp(f)/lognormal_mean;
                    fluct_coef_sum += fluct; pts++; stdev += std::pow(1.0-fluct,2.0);
                }
                else
                    fluct=0;
                row.push_back(fluct);
                
            }
            qs_fluctuation.push_back(row);
        }
        
        cout << "# Sampled local Q_s^2 fluctuations, grid " << size << "x" << size << ", cell size " << 2.0*size/(points-1) / FMGEV << " fm, fluctuation width " << Qs_fluctuation_sigma << ", mean of lognormal distribution " << lognormal_mean <<", average Q_s^2 modification " << fluct_coef_sum/pts << ", average dev: sqrt(<(modification-1)^2>) = " << std::sqrt(stdev/pts) << endl;
    }
    else
    {
        cerr << "Unknown fluctuation type set!" << endl;
        exit(1);
    }
    
}

/*
 * Get Q_s fluctuation of the given quark
 */
double Ipsat_Proton::GetQuarkQsFluctuation(unsigned int i)
{
    if (i >= qs_fluctuations_quarks.size())
    {
        cout << "Warning, asked Q_s fluctuations for quark " << i << "/" <<qs_fluctuations_quarks.size() << endl;
        return 1.0;
    }
    return qs_fluctuations_quarks[i];
}


/*
 * Get Q_s fluctuations on the transverse plane
 * Note that if we have Q_s fluctuations individually for each quark, then these
 * fluctuations are taken care in QuarkThickness() function
 */
double Ipsat_Proton::GetQsFluctuation(double x, double y)
{
    if (fluctuation_shape != LOCAL_FLUCTUATIONS or std::abs(Qs_fluctuation_sigma)<1e-5 or qs_fluctuation_coordinates.size()<1)
        return 1.0;
    
    int xind;
    if (x < qs_fluctuation_coordinates[0] or x > qs_fluctuation_coordinates[qs_fluctuation_coordinates.size()-1])
        return 1.0;
    else
        xind = FindIndex(x, qs_fluctuation_coordinates);
    int yind;
    if (y < qs_fluctuation_coordinates[0] or y > qs_fluctuation_coordinates[qs_fluctuation_coordinates.size()-1])
        return 1.0;
    else
        yind = FindIndex(y, qs_fluctuation_coordinates);
    
    //cout << x << " " << y << " " <<std::exp(qs_fluctuation[yind][xind]) << endl;
    
    return qs_fluctuation[yind][xind];
    
}

void Ipsat_Proton::SetQsFluctuation(double s)
{
    Qs_fluctuation_sigma = s;
}

void Ipsat_Proton::Init()
{
    B_p = 0.0; // GeV^-2
    B_q = 4.0;
    shape = GAUSSIAN;
    skewedness=false;
    Qs_fluctuation_sigma=0;
    fluctuation_shape = FLUCTUATE_QUARKS;
    proton_structure = QUARKS ;
    fluxtube_normalization = -1;
    origin_at_center_of_mass = false;
    number_of_quarks=3;
    
    intworkspace_zint = gsl_integration_workspace_alloc(INTPOINTS_ZINT);
}
Ipsat_Proton::Ipsat_Proton()
{
    
    allocated_gdist = false;
    
#ifdef USE_FORTRAN_IPSAT12
    ipsat = IPSAT12;
#else
    ipsat = MZSAT;
#endif
    saturation=true;
    
                                           // (4.939286653112, 1.1, 0.009631194037871, 3.058791613883, 1.342035015621);
    Init();
    
}
Ipsat_Proton::Ipsat_Proton(Ipsat_version version)
{
    allocated_gdist = false;
    
    if (version == MZSAT)
    {
        ipsat = MZSAT;
        double C=2.2894; double mu0 = std::sqrt(1.1); double lambdag=0.08289; double Ag=2.1953; double mc=1.3528;
        mzipsat = new MZ_ipsat::DipoleAmplitude(C, mu0, lambdag , Ag , mc );
        mzipsat->SetSaturation(true);
        saturation = true;
    }
    else if (version == MZNONSAT)
    {
        double C = 4.2974; double mu0=std::sqrt(1.1); double lambdag = -0.006657; double Ag=3.0391; double mc = 1.3504;
        ipsat=MZNONSAT;
        mzipsat = new MZ_ipsat::DipoleAmplitude(C, mu0, lambdag, Ag, mc);
        mzipsat->SetSaturation(false);
        saturation=false;
        
    }
    else if (version == IPSAT12)
    {
#ifndef USE_FORTRAN_IPSAT12
        cerr << "IPsat12 support is not complied (requires Fortran compiler!)" << endl;
        exit(1);
#endif
        ipsat = IPSAT12;
        saturation=true;
    }
    else if (version == IPSAT06)
    {
        gdist = new DGLAPDist();
        allocated_gdist=true;
        saturation = true;
        
    }
#ifdef USE_LCPT_DIPOLE
    else if (version == LCPT)
    {
        lcpt_dipole = new LCPT_Dipole("/Users/hejajama/Nextcloud/projects/rhorho/dipole_2d_data/x_0.01/fixed_ir_nlo_mc_5e7_mq_0.2_as_0.25_large.dat");
        lcpt_dipole->Set_out_of_range_warnings(false);
        ipsat = LCPT;
        saturation=true;
    }
#endif
    Init();
    
}
Ipsat_Proton::Ipsat_Proton(DGLAPDist *gd)
{
    gdist = gd;
    allocated_gdist = false;
    ipsat = IPSAT06;
    saturation=true;

    Init();
}

Ipsat_Proton::~Ipsat_Proton()
{
    if (allocated_gdist)
        delete gdist;
    
    if (ipsat == MZSAT or ipsat==MZNONSAT)
        delete mzipsat;
    gsl_integration_workspace_free(intworkspace_zint);
    
#ifdef USE_LCPT_DIPOLE
    if (ipsat == LCPT)
        delete lcpt_dipole;
#endif
}




/* extract collinear factorization gluon distribution from ipsat
 */
double Ipsat_Proton::xg(double x, double r)
{
    if (ipsat == IPSAT06)
    {
        double gd = gdist->Gluedist(x, r*r);
        double musqr = 4.0/(r*r)+1.17;
        double as = 12.0*M_PI/( (33.0-2.0*3.0)*log(musqr/(0.2*0.2) ) );
        
        return 2.0*3.0/(M_PI*M_PI*as ) * gd;
    }
    else if (ipsat == IPSAT12)
    {
#ifndef USE_FORTRAN_IPSAT12
        cerr << "IPsat12 support is not complied (requires Fortran compiler!)" << endl;
        exit(1);
#else
        double tmpb = 0;
        
        double n = dipole_amplitude_(&x, &r, &tmpb, &IPSAT12_PAR)/2.0;
        double exp = std::log(1.0-n);
        // exp is -pi^2 r^2/(2Nc) as xg T(0)
        double musqr = 4.0/(r*r) + 1.428;   // 1.428 corresponds to m_c=1.4 GeV
        if (IPSAT12_PAR != 2)
        {
            cerr << "Warning: mc=1.4 parameters used when extracting xg, but we have mc=1.27 IPsat! " << LINEINFO << endl;
        }
        double as = 12.0*M_PI / ( (33.0-2.0*3.0)*log(musqr/(0.156*0.156) ) );
        
        return -2.0*3.0/ (M_PI*M_PI * as* r*r * 1.0/(2.0*M_PI*4.0))*exp;
#endif
    }
    else if (ipsat == MZSAT or ipsat==MZNONSAT)
    {
        return mzipsat->xg(x, mzipsat->MuSqr(r));
    }
    else
    {
        cerr << "Ipsat version not implemented ffor Ipsat_Proton::xg" << endl;
    }
    
}

/*
 * d ln xg / d ln (1/x)
 */
double dhelperf_xg(double y, void* p);
struct dhelper_xg { Ipsat_Proton* proton; double r; };
double Ipsat_Proton::LogDerivative_xg(double x, double r)
{
    gsl_function F;
    F.function=&dhelperf_xg;
    dhelper_xg par; par.proton=this; par.r=r;
    F.params = &par;
    double result,abserr;
    double y = std::log(1.0/x);
    gsl_deriv_central (&F, y, 0.1 , &result, &abserr);
    
    return result;

}

double dhelperf_xg(double y, void* p)
{
    dhelper_xg *par = (dhelper_xg*)p;
    double x = std::exp(-y);
    return std::log(par->proton->xg(x, par->r));
}


std::vector<Vec> &Ipsat_Proton::GetQuarks()
{
    return quarks;
}

std::vector<double> Ipsat_Proton::GetRadii()
{
    std::vector<double> radii;
    for (int i=0; i<quark_bp.size(); i++)
        radii.push_back(std::sqrt(2.0*quark_bp[i]));
    return radii;
}

// Helpers used to integrate z direction of the exponential density profile
struct inthelper_exponential3d
{
    Ipsat_Proton* proton;
    double r;
    double bq;
};
double inthelperf_exponential3d(double z, void *p)
{
    inthelper_exponential3d* par = (inthelper_exponential3d*)p;
    return 1.0/(8.0*M_PI*std::pow(par->bq,3.0)) * std::exp( - std::sqrt(par->r*par->r + z*z) / par->bq);
}

double Ipsat_Proton::QuarkThickness(double r, int i)
{
    if (shape == GAUSSIAN or shape == ALBACETE or shape==MORELAND)
    {
        double bp = quark_bp[i];
        double fluct = 1.0;
        if (fluctuation_shape == FLUCTUATE_QUARKS)
        {
            fluct = GetQuarkQsFluctuation(i);
        }
        return fluct/(2.0*M_PI*bp)*std::exp(- r*r / (2.0*bp));
    }
    else if (shape == EXPONENTIAL)
    {
        double bp = quark_bp[i];
        inthelper_exponential3d par; par.proton=this; par.r=r; par.bq = bp;
        gsl_function f; f.function=&inthelperf_exponential3d;
        f.params=&par;
        double result,error;
        //gsl_integration_workspace * w =  gsl_integration_workspace_alloc (7);
        gsl_integration_qag (&f, 0, 999, 0, 5e-2, INTPOINTS_ZINT, GSL_INTEG_GAUSS15,
                             intworkspace_zint, &result, &error);
        //gsl_integration_workspace_free(w);

        
        double fluct = 1.0;
        if (fluctuation_shape == FLUCTUATE_QUARKS)
        {
            fluct = qs_fluctuations_quarks[i];
        }
        
        return 2.0*result * fluct;  // Factor 2 as we integrate from 0 to inf
        
        // 2d exponential
        //return fluct/(2.0*M_PI*bp*bp)*std::exp(- r / bp);
    }
    else
    {
        cerr << "Constituent shape unknown! " << LINEINFO << endl;
        exit(1);
    }
}

/*
 * Quark distances from the origin are sampled from this distribution
 */
double Ipsat_Proton::GaussianRadiusDistribution(double r)
{
    return std::sqrt(std::exp(1) / std::exp(B_p))*r*std::exp( - r*r / (2.0*B_p));
}


double Ipsat_Proton::ExponentialDistribution(double x, double y, double z)
{
    //a = \sqrt{12}/R_p = 3.87, with R_p = 0.895 from for 3d exp(-a*r)
    //http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1
    
    return std::exp( - std::sqrt( x*x + y*y + z*z ) / B_p );
    
}
/*
 * Distribution from 1605.09176
 *  prod i<j 1..3 (1-exp( -|s_i-s_j|^2/rc^2))
 */
double Ipsat_Proton::RepulsiveGaussianDistribution(std::vector<Vec> quarks, double rc)
{
    double result = 1.0;
    for (unsigned int i=0; i<quarks.size(); i++)
    {
        for (unsigned int j=i+1; j<quarks.size(); j++)
        {
            Vec d = quarks[i] - quarks[j];
            result *= 1.0 - std::exp( - d.LenSqr() / (rc*rc)) ;
        }
    }
    
    return result;
    
}
struct inthelper_fluxtube_z{ Ipsat_Proton* proton; Vec b;  std::vector<Vec> quarks; Vec center;};
double inthelperf_fluxtube_z(double z, void* p);

double Ipsat_Proton::FluxTubeThickness(Vec b)
{
    // Connect quarks with flux tubes that merge at the center
    // Idea from hep-lat/0606016
    // Start tubes from the quarks, and they merge at the Fermat point of the triangle
    // formed by the quarks
    // In the Ref. that point is actually Fermat point, but make this simpler now
    // The tube density is then assumed to be a Gaussian function of the distance from a line
    // connecting the center and a quark. As there are in general 3 distances, the smallest distance is chosen
    
    // First we project the tubes to z=0 and find which tube is closest. Then, integrate
    // that tube over z coordinate
    
    // Check that we are normalized
    if (fluxtube_normalization < 0)
        NormalizeFluxTubeThickness();
    

    double dist = 99999999;
    int min_index=-1;
    // Calculate distances from every line
    
    // Code to calculate projected density
    /*
    for (unsigned int i=0; i<quarks.size(); i++)
    {
        Vec line = center-quarks[i];
        Vec point = b - quarks[i];  // From quark to point b = origin moved to quark i
        
        // If point is further away from the quark or the center point than the distance between the quark and the center, then
        // we are at the edge of the tube and make it decay as a Gaussian
        // First calculate distance to the center
        Vec point_c = center - b;
        // Check
        if (point.LenSqr() > line.LenSqr() or point_c.LenSqr() > line.LenSqr())
        {
            double mind = std::min(point.Len(), point_c.Len());
            
            if (mind < dist)
            {
                dist=mind;
                min_index = i;
            }
            continue;
        }

        double proj = line*point;
        double norm = line*line;
        Vec projvec = line;
        projvec*=proj/norm;
        
        Vec distvec = point - projvec;
        
        //cout << "Distance of point " << b << " from line " << i << " is " << distvec.Len() << endl;
        
        if (distvec.Len() < dist)
        {
            dist = distvec.Len();
            min_index=i;
        }
    }*/
    
    // Closest tube is between quarks[min_index] and center.
    // That we integrate over z at point b
    inthelper_fluxtube_z par; par.b=b;
    par.proton=this;
    par.center = center3d;
    par.quarks = quarks3d;
    gsl_function f; f.params=&par;
    f.function = &inthelperf_fluxtube_z;
    double result,error;
    //gsl_integration_workspace * w =  gsl_integration_workspace_alloc (10);
    gsl_integration_qag (&f, -99, 99, 0, 1e-2, INTPOINTS_ZINT, GSL_INTEG_GAUSS15,
                          intworkspace_zint, &result, &error);
    //gsl_integration_workspace_free(w);
    
    //cout << QuarkThickness(dist, 0) << " " << result;
    
    return fluxtube_normalization*result;
    
    //return fluxtube_normalization * QuarkThickness(dist, 0);
    //cout << y<< " " << x << " " << tpsum << endl;
}

double inthelperf_fluxtube_z(double z, void* p)
{
    inthelper_fluxtube_z* par = (inthelper_fluxtube_z*)p;
    
    Vec b3d(par->b.GetX(), par->b.GetY(), z);
    
    /*cout << "Quarks: " << endl;
    cout << par->quarks[0] << endl << par->quarks[1] << endl << par->quarks[2] << endl;
    cout << "center " << par->center << endl;
    cout << "b " << b3d << endl;
    */
    double mindist=999999999999;
    for (unsigned int i=0; i<par->quarks.size(); i++)
    {
    
        // Move origin to center
        Vec quark_to_b = b3d - par->quarks[i];
        Vec quark_to_center = par->center - par->quarks[i];
        Vec center_to_b = b3d - par->center;
        
        
        
        // If we are "outside" the line, then density decreases as a Gaussian
        // from the quark/center
        // We also end up here if one angle of the triangle is larger than 120 degrees, when the Fermat
        // point is actually at one of the quarks

        if (quark_to_b.LenSqr() > quark_to_center.LenSqr() or center_to_b.LenSqr() > quark_to_center.LenSqr())
        {
            //cout << "distance from quark "<< i << ": " << quark_to_b.Len() << "   quark_to_center dist " << quark_to_center.Len() << endl;
            double dist = std::min( quark_to_b.Len(), center_to_b.Len());
            if (dist < mindist)
            {
                mindist=dist;
            }
            continue;
            //return par->proton->QuarkThickness(dist, 0);
        }
        
        // Calculate distance from the tube
        double projection_dotprod = quark_to_b * quark_to_center;
        double scaling = projection_dotprod / quark_to_center.LenSqr();
        Vec projection = quark_to_center;
        projection*= scaling;
        
        Vec dist = quark_to_b - projection;
        

        if (dist.Len() < mindist)
            mindist = dist.Len();
        
    }
    
    //cout << "mindist " << mindist << " b " << par->b << endl;

    return par->proton->QuarkThickness(mindist, 0);
    
    
}

struct inthelper_fluxtube { Ipsat_Proton* proton; double y; };
double inthelperf_fluxtube_y(double y, void* p);
void Ipsat_Proton::NormalizeFluxTubeThickness()
{
    //cout << "Normalizing fluxtube\n" << endl;
    //fluxtube_normalization = 1.0;
    //return;
    
    /// NOTE: not used, as we acutally want total energy to depend on normalization
    // FluxTubeThicknes should be normalized to unity, so calculate
    // \int dx dy FluxTubeThickness(x,y)
    // set normalization factor to 1 for this calculation
    
    
    fluxtube_normalization = 1.0;
    gsl_function f;
    inthelper_fluxtube par; par.proton = this;
    f.params = &par;
    f.function = &inthelperf_fluxtube_y;
    double result,error;
    gsl_integration_workspace * w =  gsl_integration_workspace_alloc (10);
    gsl_integration_qag (&f, -99, 99, 0, 1e-2, 10, GSL_INTEG_GAUSS15,
                          w, &result, &error);
    gsl_integration_workspace_free(w);
    //cout << "Fluxtube normalization " << result << " pm " << error << endl;
    
    cout << result << endl;;
    exit(1);
    //fluxtube_normalization = 1.0/result;
    
    
    
}
double inthelperf_fluxtube_x(double x, void* p);
double inthelperf_fluxtube_y(double y, void* p)
{
    inthelper_fluxtube* par = (inthelper_fluxtube*)p;
    par->y = y;
    gsl_function f;
    f.function = inthelperf_fluxtube_x;
    f.params = par;
    double result,error;
    gsl_integration_workspace * w =  gsl_integration_workspace_alloc (10);
    gsl_integration_qag (&f, -99, 99, 0, 1e-2, 10, GSL_INTEG_GAUSS15,
                          w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}
double inthelperf_fluxtube_x(double x, void* p)
{
    inthelper_fluxtube* par = (inthelper_fluxtube*)p;
    Vec b(x, par->y);
    return par->proton->FluxTubeThickness(b);
}

std::string Ipsat_Proton::InfoStr()
{
    std::stringstream ss;
    ss << "#IPsat proton consists of quarks at coordinates " << endl;;
    for (int i=0; i<quarks.size(); i++)
    {
        ss << "# (" << quarks3d[i].GetX() << ", " << quarks3d[i].GetY() << ", " << quarks3d[i].GetZ() << "), r=" << std::sqrt(2.0*quark_bp[i]) <<", B_q=" << quark_bp[i]<< endl ;
    }
    ss << "# Proton (Gaussian) radius " << std::sqrt(2.0*B_p) << " GeV^-1, B_p=" << B_p << " ";
    if (shape == GAUSSIAN)
        ss << "Gaussian distribution exp(-b^2/(2B_p))";
    else if (shape == EXPONENTIAL)
        ss << "Exponential distribution, exp(-b/B)";
    else if (shape==MORELAND)
        ss << "Scott Moreland's protons";
    
    ss << endl;
    
    ss << "# Structure: ";
    if (proton_structure == QUARKS)
        ss << "quarks";
    else if (proton_structure == CENTER_TUBES)
        ss << "color tubes merging at center, normalization " << fluxtube_normalization;
    if (origin_at_center_of_mass)
        ss << "Origin moved to the center of mass (shrinks proton size!)";
    
    ss << endl << "# Saturation: ";
    if (saturation)
        ss << " enabled";
    else ss << "disabled";
    ss << endl <<"# ";
    if (ipsat == IPSAT06)
        ss << "IPsat version: 2006 (arXiv:hep-ph/0304189)" << endl;
    else if (ipsat == IPSAT12)
        ss << "IPsat version: 2012 (arXiv:1212.2974)" << endl;
    else if (ipsat == MZSAT or ipsat==MZNONSAT)
        ss << "MZipsat fit (arXiv:1804.05311)" << endl;
    else if (ipsat == LCPT)
        ss << "LCPT by Dumitru, Mäntysaari, Paatelainen arXiv:2103.11682 [hep-ph]" << endl;;
    ss << "# Skewedness in dipole amplitude: ";
    if (skewedness)
        ss << " Enabled";
    else
        ss << "Disabled";
    ss << endl << "# ln Q_s^2 fluctuation width: " << Qs_fluctuation_sigma << endl;
    return ss.str();
}


void Ipsat_Proton::SetProtonWidth(double bp)
{
    B_p = bp;
}

void Ipsat_Proton::SetQuarkWidth(double bq)
{
    B_q = bq;
}

void Ipsat_Proton::SetShape(Proton_shape s)
{
    shape = s;
}

void Ipsat_Proton::SetFluctuationShape(Fluctuation_shape s)
{
    fluctuation_shape = s;
}

Fluctuation_shape Ipsat_Proton::GetFluctuationShape()
{
    return fluctuation_shape;
}

double Ipsat_Proton::Amplitude(double xpom, Vec q1, Vec q2)
{
    double quark[2] = {q1.GetX(), q1.GetY() };
    double antiquark[2] = {q2.GetX(), q2.GetY() };
    return Amplitude(xpom, quark, antiquark);
}

double Ipsat_Proton::Skewedness(double lambda)
{
    if (lambda+4.0 > GSL_SF_GAMMA_XMAX)
    {
        std::cerr << "Cant calculate skewedness, overflow, lambda=" << lambda << std::endl;
        return 1.0;
    }
    gsl_sf_result gamma_lambda_52;
    int status1 =gsl_sf_gamma_e(lambda+5.0/2.0, &gamma_lambda_52);
    gsl_sf_result gamma_lambda_4;
    int status2 =gsl_sf_gamma_e(lambda+4.0, &gamma_lambda_4);
    if (status1 or status2)
    {
        std::cerr << "Gamma function evaluation failed, Gamma(" <<lambda+5.0/2.0 << ")=" << gamma_lambda_52.val <<", Gamma(" << lambda+4.0 << ")=" << gamma_lambda_4.val << ", lambda=" << lambda << std::endl;
        return 1.0;
    }
    return std::pow(2.0, 2.0*lambda+3.0)/std::sqrt(M_PI) * gamma_lambda_52.val/gamma_lambda_4.val;
}

void Ipsat_Proton::SetSkewedness(bool s)
{
    skewedness = s;
}

void Ipsat_Proton::SetStructure(Structure s)
{
    proton_structure = s;
}

void Ipsat_Proton::SetFluxTubeNormalization(double n)
{
    fluxtube_normalization = n;
}

void Ipsat_Proton::SetQuarkCenterOfMassToOrigin(bool s)
{
    origin_at_center_of_mass = s;
}

double Ipsat_Proton::Amplitude_bint(double xpom, double r)
{
    if (ipsat == MZSAT or ipsat==MZNONSAT)
        return mzipsat->N_bint(r, xpom);
    
    cerr << "Amplitude_bint only impelmented for MZfit (ipsat_version must be " << MZSAT << " or " << MZNONSAT << ", I got " << ipsat << ")" << endl;
    return 0;
        
}

double Ipsat_Proton::Amplitude_sqr_bint(double xpom, double r)
{
    if (ipsat == MZSAT or ipsat==MZNONSAT)
        return mzipsat->N_sqr_bint(r, xpom);
    
    cerr << "Amplitude_bint_sqr only impelmented for MZfit" << endl;
    return 0;
    
}


