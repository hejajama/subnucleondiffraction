// Virtual class that can be use to easily switch between IPsat and IPGlasma
#include "src/dipole.hpp" 
// IP-Glasma
#include "src/ipglasma.hpp"
// IPsat
#include "src/ipsat_proton.hpp"

// IPsat_proton needs global random number generator
#include <gsl/gsl_rng.h>
gsl_rng* global_rng;

// 2d vectors
#include "src/vector.hpp"

using namespace std;

int main()
{
	// Initialize global random number generator
    const gsl_rng_type * rngtype = gsl_rng_default;
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(rngtype);


	bool use_ipglasma = true;	// if false, use IPsat
	DipoleAmplitude* amplitude;	// This will be either IP-Glasma or IPsat

	// If we use IP-sat
	if (use_ipglasma == false)
    {
		amplitude = new Ipsat_Proton;
		((Ipsat_Proton*)amplitude)->SetProtonWidth(0);
		((Ipsat_Proton*)amplitude)->SetQuarkWidth(4);	// no fluctuations	
		// Proton "size" is B_p=4.0 Gev^(-2), without fluctuations we set three quarks/hot spots
		// at the center (proton width 0), and set the width of the quarks to be 4
		// This is the IPsat parametrization used to fit the HERA data by Amir et al
	}
	else 
	{
		// Load wilson line, the 2nd argument is the step size of the grid
		// Note that this is the only place where a dimensionful number is in fm, not GeV
		// So this means that in the IP-Glasma grid lattice spacing a=0.01 fm.
		amplitude = new IPGlasma("wilsonline.txt", 0.01);
		
		// If you have slightly older version of the code, the IPGlasma constructor
		// does not support the step size parameter, but the step is set in the IPGlasma::IPGlasma()
		// methdod.	But it is easy to update the local copy of the git repository
	}

	amplitude->InitializeTarget();

	double xbj=0.01;   // IPsat supports Bjorken-x, in IP-Glasma this is
						// 0.01 * exp(-rapidity from input file) 
	Vec impact_param(5.068, 0);	// Note: all dimensionful numbers are GeV^n

	// dipole size
	Vec r(-1,2);	
	
	// Positions for q and antiq
	Vec q1 = r*0.5 +impact_param;
	Vec q2 = r*(-0.5) + impact_param;
	cout << "Dipole amplitude, real part, quar positions:" << endl; 
	cout << q1 << endl;
	cout << q2 << endl;
	cout << amplitude->Amplitude(xbj, q1, q2) << endl;
	cout << "And imaginary part:" << endl;
	cout << amplitude->AmplitudeImaginaryPart(xbj, q1, q2) << endl;


	delete amplitude;
}
	
