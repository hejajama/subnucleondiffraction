/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a IPglasma nucleus
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */


#include "ipglasma.hpp"
#include "subnucleon_config.hpp"
#include <string>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <sstream>
#include <cstdlib>

typedef unsigned int uint;

using std::cout;
using std::cerr;
using std::endl;
extern gsl_rng *global_rng;

const int NC=3;

int FindIndex(double val, std::vector<double> &vec);


/*
 * Calculate dipole amplitude form Wilson lines
 * Search the closest grid point that corresponds to the given quark/antiquark coordinate
 * Then calcualte 1 - 1/Nc Tr U(quark) U^dagger(antiquark)
 */
double IPGlasma::Amplitude(double xpom, double q1[2], double q2[2] )
{
    ApplyPeriodicBoundaryConditions(q1);
    ApplyPeriodicBoundaryConditions(q2);
    
    // Out of grid? Return 0 (probably very large dipole)

    if (q1[0] < xcoords[0] or q1[0] > xcoords[xcoords.size()-1]
        or q1[1] < ycoords[0] or q1[1] > ycoords[ycoords.size()-1]
        or q2[0] < xcoords[0] or q2[0] > xcoords[xcoords.size()-1]
        or q2[1] < ycoords[0] or q2[1] > ycoords[ycoords.size()-1])
    {
        if (periodic_boundary_conditions)
            cerr << "WTF, I'm here..." << endl;
            
        return 0;
        
    } 

    double  r = sqrt( pow(q1[0]-q2[0],2) + pow(q1[1]-q2[1],2));

	// Smaller dipole than the grid
	if (r < std::abs( xcoords[1] - xcoords[0])) 
			return 0;	
	
    if (schwinger and r > schwinger_rc)
    {
	    if (q1[0] < xcoords[0]) q1[0] = xcoords[0];
	    if (q1[1] < xcoords[0]) q1[1] = xcoords[0];
	    if (q2[0] < xcoords[0]) q2[0] = xcoords[0];
	    if (q2[1] < xcoords[0]) q2[1] = xcoords[0];
	    if (q1[0] > xcoords[xcoords.size()-1]) q1[0] = xcoords[xcoords.size()-1];
	    if (q2[0] > xcoords[xcoords.size()-1]) q2[0] = xcoords[xcoords.size()-1];
	    if (q1[1] > xcoords[xcoords.size()-1]) q1[1] = xcoords[xcoords.size()-1];
            if (q2[1] > xcoords[xcoords.size()-1]) q2[1] = xcoords[xcoords.size()-1]; 

            double b = sqrt( pow( (q1[0]+q2[0])/2.0, 2.0) + pow((q1[1]+q2[1])/2,2.0) );
	    
            // Dipole schwinger
	   
            // Kind of schwinger
 	    double q3[2];
	
	    // Put new quark randomly between q1 and q2, calculate weighted
	    // mean with weign in [0.3, 0.7]
	    double w = gsl_rng_uniform(global_rng)*0.7+0.3;

	    q3[0] = (w*q1[0]+(1.0-w)*q2[0])/1.0;
	    q3[1] = (w*q1[1]+(1.0-w)*q2[1])/1.0;

	    double s1 = 1.0 - Amplitude(xpom, q1, q3);
	    double s2 = 1.0 - Amplitude(xpom, q2, q3);
	    if (s1 < 0) s1 = 0; if (s1>1) s1=1;
	    if (s2<0) s2=0; if (s2>1) s2=1;
	    //cout << "r " << r << " Schwinger gives " << s1 << " and " << s2 << endl;
	    return 1.0 - s1*s2;
		
	  // Quadrupole schwinger
/* 
	   int splits = 2 * int(r / schwinger_rc) + 1;
           double w = 1.0 / ((double)splits);
		cout << "test, rpoints (" << q1[0] << ", " << q1[1] << ") and (" << q2[0] << ", " << q2[1] << ")" << endl;
	   std::vector<WilsonLine> wlines;
           for (unsigned int i=0; i<splits+1; i++)
           {
		double xy[2];
		xy[0] = i*w*q2[0] + (1.0 - i*w)*q1[0];
		xy[1] = i*w*q2[1] + (1.0 - i*w)*q1[1];
		cout << "New point " << xy[0] << " " << xy[1] << endl;
		WilsonLine tmp = GetWilsonLine(xy[0], xy[1]);
		if (i % 2 == 1)
			tmp = tmp.HermitianConjugate();
		wlines.push_back(tmp);
	    }
            WilsonLine r;
	    r = wlines[0];
	    for (int i=1; i<wlines.size(); i++)
{
	//cout << "Product " << i << endl;
		r = r * wlines[i];	
}	

	cout << "Trace " <<  1.0 - r.Trace().real()/NC << endl;
*/


/*
	  
	    double q3[2]; double q4[2];
	    q3[0] = 0.333333 * q2[0] + (1.0 - 0.333333) * q1[0];
	    q3[1] = 0.333333 * q2[1] + (1.0 - 0.333333) * q1[1];
            q4[0] = 0.666666 * q2[0] + (1.0 - 0.666666) * q1[0];
	    q4[1] = 0.6666666 * q2[1] + (1.0 - 0.666666) * q1[1];

	    WilsonLine w1 = GetWilsonLine(q1[0], q1[1]);
            WilsonLine w2 = GetWilsonLine(q2[0], q2[1]);
	    w2 = w2.HermitianConjugate();
            WilsonLine w3 = GetWilsonLine(q3[0], q3[1]);
            WilsonLine w4 = GetWilsonLine(q4[0], q4[1]);
	    w3 = w3.HermitianConjugate();

	    WilsonLine p;
            p = w1 * w3 * w4 * w2;

	//	cout << "old points (" << q1[0] << ", " << q1[1] << "), (" << q3[0] << ", " << q3[1] << "), and (" << q4[0] << ", " << q4[1] << ")  (" << q2[0] << ", " << q2[1] << ")";
	//	cout << " trace " << 1.0 - p.Trace().real()/NC << endl;
            return 1.0 - p.Trace().real()/NC;

	*/		
		
     }
    // First find corresponding grid indeces
    WilsonLine quark = GetWilsonLine(q1[0], q1[1]);
    WilsonLine antiquark = GetWilsonLine(q2[0], q2[1]);
    
    //antiquark = antiquark.HermitianConjugate();
    
    WilsonLine prod;
    try {
        prod = quark.MultiplyByHermitianConjugate(antiquark);
        //prod =  quark*antiquark;
    } catch (...) {
        cerr << "Matrix multiplication failed!" << endl;
        cout << "Quark: " << q1[0] << ", " << q1[1] << endl;
        cout << quark << endl;
        cout << "Antiquark: " << q2[0] << ", " << q2[1] << endl;
        cout << antiquark << endl;
        exit(1);
    }
    std::complex<double > amp =  1.0 - 1.0/NC * prod.Trace();
    
    double result = amp.real();
    if (result < 0) return 0;
    return result;
    if (result > 1)
        return 1;
    if (result < 0)
        return 0;

    if (isnan(result))
    {
        cerr << "Wilson line trance NaN, quark coords " << q1[0] << ", " << q1[1] << " and " << q2[0] << ", " << q2[1] << endl;
	exit(1);
    } 
     
    return result;
}
// Stupid copypaste
double IPGlasma::AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] )
{

    ApplyPeriodicBoundaryConditions(q1);
    ApplyPeriodicBoundaryConditions(q2);

    // Out of grid? Return 1 (probably very large dipole)
    if (q1[0] < xcoords[0] or q1[0] > xcoords[xcoords.size()-1]
        or q1[1] < ycoords[0] or q1[1] > ycoords[ycoords.size()-1]
        or q2[0] < xcoords[0] or q2[0] > xcoords[xcoords.size()-1]
        or q2[1] < ycoords[0] or q2[1] > ycoords[ycoords.size()-1])
        return 0;
double  r = sqrt( pow(q1[0]-q2[0],2) + pow(q1[1]-q2[1],2));
        if (r < std::abs( xcoords[1] - xcoords[0]))
                        return 0;
    
    // First find corresponding grid indeces
    WilsonLine quark = GetWilsonLine(q1[0], q1[1]);
    WilsonLine antiquark = GetWilsonLine(q2[0], q2[1]);
//    antiquark = antiquark.HermitianConjugate();
    
    WilsonLine prod;
    try {
        //prod =  quark*antiquark;
        prod = quark.MultiplyByHermitianConjugate(antiquark);

    } catch (...) {
        cerr << "Matrix multiplication failed!" << endl;
        cout << "Quark: " << q1[0] << ", " << q1[1] << endl;
        cout << quark << endl;
        cout << "Antiquark: " << q2[0] << ", " << q2[1] << endl;
        cout << antiquark << endl;
        exit(1);
    }
    std::complex<double > amp =  1.0 - 1.0/NC * prod.Trace();
    
    double result = amp.imag();
    return result;
    if (result > 1)
        return 1;
    if (result < 0)
        return 0;
    
    return result;
}

WilsonLine& IPGlasma::GetWilsonLine(double x, double y)
{
    double q[2]={x,y};
    ApplyPeriodicBoundaryConditions(q);
    x=q[0];
    y=q[1];

    std::vector<int> coords = LatticeCoordinates(x,y);

    // Handle edges
    if (coords[0] < 0)
        coords[0] = 0;
    if (coords[0] >= xcoords.size())
        coords[0] = xcoords.size()-1;

    if (coords[1] < 0)
        coords[1] = 0;
    if (coords[1] >= ycoords.size())
        coords[1] = ycoords.size()-1;
    
    return GetWilsonLine( WilsonLineCoordinate(coords[0],coords[1]));
    
}

IPGlasma::IPGlasma(std::string file)
{
    int load = LoadData(file, 5.12/512.0);     // By default assume lattice spacing 0.01fm, or 5.12fm 512^2 lattice

    if (load<0)
        exit(1);

    periodic_boundary_conditions=false;
}

IPGlasma::IPGlasma(std::string file, double step, WilsonLineDataFileType type)
{
   int load =LoadData(file, step, type);
    
    if (load<0)
        exit(1);

    periodic_boundary_conditions=false;
}


// Note that here step is in fm!
int IPGlasma::LoadData(std::string fname, double step, WilsonLineDataFileType type)
{
    // Load data
    datafile=fname;
    
    if (type == BINARY)
        return LoadBinaryData(fname, step);
    
    // Syntax: x y [fm/indeces] matrix elements Re Im for elements (0,0), (0,1), (0,2), (1,0), ...
    std::ifstream f(fname.c_str());
    
    if (!f.is_open())
    {
        std::cerr << "Could not open file " << fname << " " << LINEINFO << std::endl;;
        return -1;
    }
    std::string line;
    
    while(!f.eof() )
    {
        std::getline(f, line);
        if (line[0]=='#' or line.length() < 10)   // Comment line, or empty
        {
            continue;
        }
        
        // Parse
        std::stringstream ss(line);
        // Get coordinates
        double x,y;
        ss >> x;
        ss >> y;
        
        // Datafile is in fm, but we want to use GeVs in this code
        // Once we have load all points, we will sift all coordinates such that 0 is at the center
        

        x = step*x*FMGEV;
        y = step*y*FMGEV;
        
        
        // Read rows and columns
        std::vector< std::vector< std::complex<double> > > matrix;
        for (int row=0; row < 3; row++)
        {
            std::vector< std::complex<double> > tmprow;
            for (int col=0; col<3; col++)
            {
                double real, imag;
                ss >> real;
                ss >> imag;
                std::complex<double> element (real, imag);
                tmprow.push_back(element);
            }
            matrix.push_back(tmprow);
        }
        
        // Save Wilson line
        WilsonLine w(matrix);
        w = w.Transpose();
        wilsonlines.push_back(w);
        
        if (w.Size() != 3)
        {
            cerr << "Matrix at coordinate " << x << " " << y << " has size " << w.Size() << endl;
            exit(1);
        }
        
        // We assume that grid is symmetric, so dont save same valeus multiple times
        // In the datafile the coordinates are increasing
        if (xcoords.size() == 0 or x > xcoords[xcoords.size()-1])
            xcoords.push_back(x);
        if (ycoords.size() == 0 or y > ycoords[ycoords.size()-1])
            ycoords.push_back(y);
        
    }
    
    // Shift coordinates s.t. 0 fm is at the center
    double center = xcoords[xcoords.size()/2];
    for (int i=0; i<xcoords.size(); i++)
    {
        xcoords[i] -= center;
        ycoords[i] -= center;
    }
    f.close();
    
    if (xcoords.size() != ycoords.size())
    {
        cerr << "xcoords.size() != ycoords.size(), probably uncomplete input data? Datafile " << fname << " -  " << LINEINFO << endl;
        return -1;
    }

    if (xcoords.size() < 10 or ycoords.size()<10)
    {
        cerr << "Grid size is " << xcoords.size() << " x " << ycoords.size() << ", this makes no sense! File " << fname << " - " << LINEINFO << endl;
        return -1;
    }
    

    
    // Now, given that we have a point (x,y), we can find the index xind such that
    // xcoords[xind] is closest to x, and similarly ind
    // Then, the corresponding Wilson line is
    // wilsonlines[ xcoords.size()*xind + yind]
    // Of course this is symmetric and we could just as well swap xind and yind

   SetSchwinger(false); 
   // std::cout <<"# Loaded " << wilsonlines.size() << " Wilson lines from file " << datafile << ", grid size " << xcoords.size() << " x " << ycoords.size() << " grid range [" << xcoords[0] << ", " << xcoords[xcoords.size()-1] << "]" << " step size " << xcoords[1]-xcoords[0] << " GeV^-1" << std::endl;

        
    return 0;
}

std::vector<int> IPGlasma::LatticeCoordinates(double x, double y)
{
    std::vector<int> ret;

    // Note: My lattice is from -L/2 to L/2, so I need to shift the coordinates
    x = x + xcoords[xcoords.size()-1];
    y = y + ycoords[ycoords.size()-1];
    double lattice_spacing = xcoords[1]-xcoords[0];

    int ix = x/lattice_spacing; 
    int iy = y/lattice_spacing; 

    return std::vector<int> {ix, iy}; 
}

int IPGlasma::WilsonLineCoordinate(int  xind, int yind)
{
    return xcoords.size()*xind + yind; 
}


/*
 * Load data from binary file
 */
int IPGlasma::LoadBinaryData(std::string fname, double step)
{
    std::ifstream InStream;
    InStream.precision(10);
    InStream.open(fname.c_str(), std::ios::in | std::ios::binary);
    int N;
    int Nc;
    double L,a;
    
    if(InStream.is_open())
    {
        // READING IN PARAMETERS -----------------------------------------------------------------//
        double temp;
        InStream.read(reinterpret_cast<char*>(&N), sizeof(int));
        InStream.read(reinterpret_cast<char*>(&Nc), sizeof(int));
        InStream.read(reinterpret_cast<char*>(&L), sizeof(double));
        InStream.read(reinterpret_cast<char*>(&a), sizeof(double));
        InStream.read(reinterpret_cast<char*>(&temp), sizeof(double));
        
        std::cout << "# BINARY Size is " << N << ", Nc " << Nc << ", length is [fm] " << L << ", a is [fm]" << a << "  (user specified " << step << ")" << std::endl;
        
        // Init Wilson lines - fill in later
        wilsonlines.clear();
        for (unsigned int i=0; i<N*N; i++)
        {
            WilsonLine w;
            wilsonlines.push_back(w);
        }
        
        
        // READING ACTUAL DATA --------------------------------------------------------------------//
        double ValueBuffer;
        int INPUT_CTR=0;
        double re,im;
        
        while( InStream.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
        {
            if(INPUT_CTR%2==0)              //this is the real part
            {
                re=ValueBuffer;
            }
            else                            // this is the imaginary part, write then to variable //
            {
                im=ValueBuffer;
                
                int TEMPINDX=((INPUT_CTR-1)/2);
                int PositionIndx = TEMPINDX / 9;
                int ix = PositionIndx / N;
                int iy = PositionIndx - N*ix;
                
                
                int MatrixIndx=TEMPINDX - PositionIndx*9;
                int j=MatrixIndx/3;
                int k=MatrixIndx-j*3;
                
                int indx = N*iy + ix;
                wilsonlines[indx].Set(j,k, std::complex<double> (re,im));
            }
            INPUT_CTR++;
        }
    }
    else
    {
        std::cerr << "ERROR COULD NOT OPEN FILE " << fname << std::endl;
    }
    
    InStream.close();
    
    // Fill x and y coordinate values in GeV^-1
    for (unsigned int i=0; i<N; i++)
    {
        xcoords.push_back(i*a*5.068);
        ycoords.push_back(i*a*5.068);
    }
    
    // Shift coordinates s.t. 0 fm is at the center
    double center = xcoords[xcoords.size()/2];
    for (int i=0; i<xcoords.size(); i++)
    {
        xcoords[i] -= center;
        ycoords[i] -= center;
    }
    
    SetSchwinger(false);
    return 0;
}

double IPGlasma::MinX()
{
    return xcoords[0];
}

double IPGlasma::MaxX()
{
    return xcoords[xcoords.size()-1];
}


double IPGlasma::XStep()
{
    return xcoords[1] - xcoords[0];
}


void IPGlasma::ApplyPeriodicBoundaryConditions(double q[2])
{
    
    if (periodic_boundary_conditions == false) return;

    double Lx = xcoords[xcoords.size()-1]-xcoords[0];
    double Ly = ycoords[ycoords.size()-1]-ycoords[0];
    while (q[0] < xcoords[0]) 
        q[0] = q[0] + Lx;

    while (q[0] > xcoords[xcoords.size()-1])
        q[0] = q[0] - Lx;

     while (q[1] < ycoords[0]) 
        q[1] = q[1] + Ly;

    while (q[1] > ycoords[ycoords.size()-1])
        q[1] = q[1] - Ly;


}

std::string IPGlasma::InfoStr()
{
    std::stringstream ss;
    ss << "# IPGlasma loaded from file " << datafile << " lattice " << xcoords.size() << "^2 range [" << xcoords[0]/5.068 << ", " << xcoords[xcoords.size()-1]/5.068 << "] fm" << endl ;
    if (periodic_boundary_conditions) ss << "# Periodic boundary conditions" << endl;
    if (schwinger) ss << "# schwinger mechanism included, rc=" << schwinger_rc << " GeV^-1" << endl;
    return ss.str();
}

std::vector<double> &IPGlasma::GetXCoordinates()
{
    return xcoords;
}

void IPGlasma::SetSchwinger(bool s, double rc)
{
    schwinger = s;
    schwinger_rc = rc;
}

