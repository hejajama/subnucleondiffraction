/*
 * Nucleus that consists of nucleons
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */

#ifndef nucleons_hpp
#define nucleons_hpp

#include "dipole.hpp"
#include <vector>
#include "ipsat_proton.hpp"
#include "vector.hpp"
#include "gdist_dglap.hpp"
#include "interpolation.hpp"

enum DeuteronStructure
{
    NUCLEONS,   // Independent p and n
    TUBE    // Connect nucleons by a tube
};

enum Nucleon
{
    Proton,
    Neutron
};

// Hulthen:
// Extended Hulthen: Phys. Rev. 151, 772
// WoodsSaxon is parameters from PHOBOS Glauber 1408.2549
enum DeuteronWaveFunctionType {
    Hulthen,
    ExtendedHulthen,
    WoodsSaxon,
    VMC, // https://www.phy.anl.gov/theory/research/density2/
};

enum NuclearDensity {
    NuclearDensity_WoodsSaxon,
    P_N_Densities
};


class Nucleons : public DipoleAmplitude
{
public:
    // Evaluate dipole ampltitude, qaurks at coordinates x1 and x2
    // Array points are x and y coordinates
    double Amplitude(double xpom, double q1[2], double q2[2] );
    
    void InitializeTarget();
    
    void SetHeId(int i);    // Set He3 id
    
    Nucleons(std::vector<DipoleAmplitude*> nucleons);
    ~Nucleons();
    double WS_unnorm(double r, Nucleon nucleon = Proton);
    
    std::string InfoStr();
    
    double DeuteronWaveFunction(double r);  // Deuteron wf, r is 3d vector
    
    double DeuteronTubeDensity(Vec b);  // Model deuteron as a Gaussian tube
    
    std::vector<DipoleAmplitude*> GetNucleons();
	std::vector<Vec> GetNucleonPositions() { return nucleon_positions; }
    
    DeuteronStructure GetDeuteronStructure() { return deuteron_structure; }
    void SetDeuteronStructure(DeuteronStructure d) { deuteron_structure = d; }
    void SetDeuteronWF(DeuteronWaveFunctionType wf) { DeuteronWF = wf; }
    DeuteronWaveFunctionType GetDeuteronWF() { return DeuteronWF; }

    void SetNuclearDensity(NuclearDensity d);
    
private:
    int A;
    int Z; // Used if P_N_DENSITIES
    std::vector<DipoleAmplitude*> nucleons; // Todo: change to support general proton/nucleon class?
    std::vector<Vec> nucleon_positions;
    
    double ws_delta;
    double ws_ra;
    
    DeuteronStructure deuteron_structure;
    DeuteronWaveFunctionType DeuteronWF;
    
    Interpolator VMC_interpolator;

    Interpolator proton_density_interpolator = Interpolator();
    Interpolator neutron_density_interpolator = Interpolator();

    NuclearDensity nuclear_density = NuclearDensity_WoodsSaxon;
    
    int he3_id; // which He3 configuration we use
    
};

#endif /* nucleons_hpp */
