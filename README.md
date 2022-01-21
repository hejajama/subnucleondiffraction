Subnucleondiffraction code

References:

 * H. Mäntysaari, B. Schenke, Phys.Rev. D94 (2016) no.3, 034042, arXiv:1607.01711
 * H. Mäntysaari, B. Schenke, Phys.Rev.Lett. 117 (2016) no.5, 052301, arXiv:1603.04349


This program calculates scattering amplitude for the process 
`gamma^* + A -> V + A`
where V is a vector meson (Upsilon, JPsi, rho, phi)

**Note** This code is constantly developed, and individual commits are not guaranteed to work properly. If you want to use this code in your project, it is probably good idea to first communicate directly with Heikki Mäntysaari <heikki.mantysaari@jyu.fi>


## Compiling

First one has to download and compile my "amplitudelib" package.
Get it from Github and compile using CMake:
    
    git clone https://github.com/hejajama/amplitudelib.git amplitudelib_v2
    cd amplitudelib_v2
    mkdir build
    cd build
    cmake ..
    make
 
Then, download this main program ("subnucleondiffraction"):

    cd ../..   (go back to the parent folder)
    git clone https://github.com/hejajama/subnucleondiffraction.git

Note that when you compile this, it is assumed that the
libraries generated above can be found from 
`../amplitudelib_v2/build/lib/`

So amplitudelib_v2 and subnucleondiffraction folders should be located
in the same directory.


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

Examples

    GSL_RNG_SEED=1 ./build/bin/subnucleondiffraction -dipole 1 ipsatproton 3.3 0.7 -real -Q2 0 -W 75 -mcintpoints 1e5
    
Calculates diffractive scattering amplitude (real part, imaginary part: use `-imag` instead of `-real`)
Q2  is the photon virtuality in GeV<sup>2</sup> and W is the center-of-mass energy. This woud use 10<sup>5</sup> MC integration points points are used in (adaptive) Monte Carlo integration. The MonteCarlo method (Vegas or MISER) can be selected using the `-mcint` flag (see `src/main.cpp`)

Using a heavy nucleus instead of proton, replace `1 -> 197` (Au) or any other A. For A>3 Woods-Saxon is used. Deuteron and 3He are handled separately.

Random seed is set by GSL_RNG_SEED enviromental variable. One **must** use the same RNG_SEED when calculating real and 
imaginary parts!

Round proton is "ipsatproton 0 4" (first number controls the Width of the Gaussian from which the hot spot locations are sampled, and the second number is the width of the hot spots.

Q_s fluctuations for constituent quarks are set as:

    -qsfluct 0.5 -qsfluctshape quarks 
    
where the first number is the width of the log-normal distribution sigma.

Wilson lines generated using the IPGlasma cdoe can be used instead of the IPsat dipole as follows: 
    
    -dipole 1 ipglasma FILENAME step

The step size in fm, and should be `L/N` (`L` is the lattice length, `N` number of lattice points)

Or, it is more efficient use Wilson lines in binary format

    -dipole 1 ipglasma_binary FILENAME

In this case there is no need to specify the step size.
 
The code outputs

    t   amplitude(T)  amplitude(L)

So amplitudes for different polarizations separately. All dimensionful units in this code are in GeV unless stated otherwise. Note that the user has to calculate real and imaginary part separately (but the imaginary part does not affect the coherent cross section).

The coherent cross section is then 

    dsigma/dt = 1/(16 pi ) <A>^2   [in 1/GeV^4]

Where \<A\> is the average of the amplitudes.

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
IPsat fit 

 * A. Rezaeian, M. Siddikov, M. Van de Klundert, R. Venugopalan Phys.Rev. D87 (2013) no.3, 034002 

The code is directly from their work, but compiled here to the library 
referred to as libColorDipole. This is not compiled by default, so one has to manually edit CMake files to compile it.

The NRQCD wave function reference is

* T. Lappi, H. Mäntysaari, J. Penttala,  Phys.Rev.D 102 (2020) 5, 054020, arXiv:2006.02830 [hep-ph]

The Boosted Gaussia and Gaus-LC wave function parametrizations are from

* H. Kowalski, L. Motyka, G. Watt, Phys.Rev.D 74 (2006) 074016, arXiv:hep-ph/0606272 

Deuteron and <sup>3</sup>He: see

* H. Mäntysaari, B. Schenke, Phys.Rev.C 101 (2020) 1, 015203, arXiv:1910.03297
