Subnucleondiffraction code

References:

 * H. Mäntysaari, B. Schenke, Phys.Rev. D94 (2016) no.3, 034042, arXiv:1607.01711
 * H. Mäntysaari, B. Schenke, Phys.Rev.Lett. 117 (2016) no.5, 052301, arXiv:1603.04349


This program calculates scattering amplitude for the process 
`gamma^* + A -> V + A`
where V is a vector meson (Upsilon, JPsi, rho, phi)

**Note** This code is constantly developed, and individual commits are not guaranteed to work properly. If you want to use this code in your project, it is probably good idea to first communicate directly with Heikki Mäntysaari <heikki.mantysaari@jyu.fi>

**Note 2** As of Nov/2024, the output format has changed, the code now always computes the real and imaginary parts and prints both results.


## Compiling
 
First clone this repository

    git clone https://github.com/hejajama/subnucleondiffraction.git
    
How to compile:

    mkdir build
    cd build
    cmake ..
    make

Code requires GSL version 2. CMake takes care of
finding the correct compiler flags.

## Usage 

How to run:
See `./build/bin/subnucleondiffraction -help`

User has to specify the dipole-target scattering amplitude. Supported dipole amplitudes are

 * IPsat with event-by-event fluctuating geometry
 * Wilson lines generated using the IPGlasma code (https://github.com/schenke/ipglasma)
 * JIMWLK-evolved Wilson lines also work, the evolution can be solved using (https://github.com/hejajama/jimwlk)

Examples

    GSL_RNG_SEED=1 ./build/bin/subnucleondiffraction -dipole 1 ipsatproton 3.3 0.7 -Q2 0 -W 75 -mcintpoints 1e5
    
Calculates diffractive scattering amplitude. Q2  is the photon virtuality in GeV<sup>2</sup> and W is the center-of-mass energy (again in GeV). This woud use 10<sup>5</sup> MC integration points in Monte Carlo integration. The MonteCarlo integration is done using the Suave routine from the Cuba package.

Using a heavy nucleus instead of proton, replace `1 -> 197` (Au) or any other A. For A>3 Woods-Saxon is used. Deuteron and 3He are handled separately.

Random seed is set by `GSL_RNG_SEED` enviromental variable. 

Round IPsat-proton is e.g. `ipsatproton 0 4` (first number controls the width of the Gaussian from which the hot spot locations are sampled, and the second number is the width of the hot spot. This code always uses three hotspots, edit `src/ipsat_proton.cpp` if necessary. If the center-of-mass should be moved to the origin, add `com` at the end:

    -dipole 1 ipsatproton 4.5 1.0 com

The magnitude of the Q_s fluctuations is set as

    -qsfluct 0.5

Wilson lines generated using the IPGlasma code can be used instead of the IPsat dipole as follows: 
    
    -dipole 1 ipglasma FILENAME step

The step size in fm, and should be `L/N` (`L` is the lattice length, `N` number of lattice points). Note: everywhere else this code uses GeV^n units.

It is more efficient use Wilson lines in binary format (generated using the `writeInitialWilsonLines 2` option in IP-Glasma)

    -dipole 1 ipglasma_binary FILENAME

In this case there is no need to specify the step size.
 
The code outputs the squared momentum transfer |t| and complex scattering amplitudes separately for the transverse and longitudinal photon, the syntax is

    t   transverse real, transverse imag, longitudinal real, longitudinal imag

 All dimensionful units in this code are in GeV unless stated otherwise, so amplitude is in 1/GeV^2 (cross section 1/GeV^4). 

The coherent diffractive cross section is 

    dsigma/dt = 1/(16 pi ) <A>^2   [in 1/GeV^4]

Where \<A\> is the average of the amplitudes. Note that the factor $1/(16\pi)$ is not included in this code!

Similalry the incoherent corss section is computed by replacing `<A>^2` by `<A^2> - <A>^2`

In order to calculate total gamma-p (or gamma-A) cross section, one can calulate F2 as follows

    GSL_RNG_SEED=1 ./diffraction -dipole 1 ipsatproton 0 4 -F2 Qsqr xbj

The output is
x   Q^2   F_2,light  F_2,charm   F_2,total   F_L,light  F_L,charm   F_L,tot 

## Additional features

Some branches in this repository include support for additional processes or other improvements. These include

* `photon_kt`: support finite photon_kt in ultra peripheral collisions
* `dvcs_vm_masses`: calculate azimuthal correlations between the outgoing lepton and the produced vector meson
* `dijet`: calculate diffractive dijet production

## Further references 
The default IPsat fit used is

 * H. Mäntysaari, P. Zurita,  Phys.Rev.D 98 (2018) 036002, arXiv:1804.05311 [hep-ph]

This program also contains codes to evaluate the dipole amplitude from the
IPsat fit (not build by default)

 * A. Rezaeian, M. Siddikov, M. Van de Klundert, R. Venugopalan Phys.Rev. D87 (2013) no.3, 034002 

The code is directly from their work, but compiled here to the library 
referred to as libColorDipole. This is not compiled by default, so one has to manually edit CMake files to compile it.

The NRQCD wave function reference is

* T. Lappi, H. Mäntysaari, J. Penttala,  Phys.Rev.D 102 (2020) 5, 054020, arXiv:2006.02830 [hep-ph]

The Boosted Gaussia and Gaus-LC wave function parametrizations are from

* H. Kowalski, L. Motyka, G. Watt, Phys.Rev.D 74 (2006) 074016, arXiv:hep-ph/0606272 

Deuteron and <sup>3</sup>He: see

* H. Mäntysaari, B. Schenke, Phys.Rev.C 101 (2020) 1, 015203, arXiv:1910.03297
