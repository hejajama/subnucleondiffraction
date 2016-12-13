/*
 * Diffraction at sub-nucleon scale
 * Dipole amplitude for a IPglasma nucleus
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2015
 */


#include "ipglasma.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <cstdlib>

using namespace Amplitude;

using std::cout;
using std::endl;

const int NC=3;

/*
 * Calculate dipole amplitude form Wilson lines
 * Search the closest grid point that corresponds to the given quark/antiquark coordinate
 * Then calcualte 1 - 1/Nc Tr U(quark) U^dagger(antiquark)
 */
double IPGlasma::Amplitude(double xpom, double q1[2], double q2[2] )
{
    
    // Out of grid? Return 0 (probably very large dipole)
    if (q1[0] < xcoords[0] or q1[0] > xcoords[xcoords.size()-1]
        or q1[1] < ycoords[0] or q1[1] > ycoords[ycoords.size()-1]
        or q2[0] < xcoords[0] or q2[0] > xcoords[xcoords.size()-1]
        or q2[1] < ycoords[0] or q2[1] > ycoords[ycoords.size()-1])
            return 0;
    
    // First find corresponding grid indeces
    WilsonLine quark = GetWilsonLine(q1[0], q1[1]);
    WilsonLine antiquark = GetWilsonLine(q2[0], q2[1]);
    antiquark = antiquark.HermitianConjugate();
    
    WilsonLine prod;
    try {
        prod =  quark*antiquark;
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
    return result;
    if (result > 1)
        return 1;
    if (result < 0)
        return 0;
     
    return result;
}
// Stupid copypaste
double IPGlasma::AmplitudeImaginaryPart(double xpom, double q1[2], double q2[2] )
{
    // Out of grid? Return 1 (probably very large dipole)
    if (q1[0] < xcoords[0] or q1[0] > xcoords[xcoords.size()-1]
        or q1[1] < ycoords[0] or q1[1] > ycoords[ycoords.size()-1]
        or q2[0] < xcoords[0] or q2[0] > xcoords[xcoords.size()-1]
        or q2[1] < ycoords[0] or q2[1] > ycoords[ycoords.size()-1])
        return 0;
    
    // First find corresponding grid indeces
    WilsonLine quark = GetWilsonLine(q1[0], q1[1]);
    WilsonLine antiquark = GetWilsonLine(q2[0], q2[1]);
    antiquark = antiquark.HermitianConjugate();
    
    WilsonLine prod;
    try {
        prod =  quark*antiquark;
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
    int xind = FindIndex(x, xcoords);
    int yind = FindIndex(y, ycoords);
    
    // Handle edges
    if (xind < 0)
        xind = 0;
    if (xind >= xcoords.size())
        xind = xcoords.size()-1;

    if (yind < 0)
        yind = 0;
    if (yind >= ycoords.size())
        yind = ycoords.size()-1;
    //cout << "Coordinates " << x << ", "  << y << " indeces " << xind << ", " << yind << endl;
    
    return wilsonlines[ xind*xcoords.size() + yind];
    
}

IPGlasma::IPGlasma(std::string file)
{
    // Load data
    datafile = file;
    // Syntax: x y [fm/indeces] matrix elements Re Im for elements (0,0), (0,1), (0,2), (1,0), ...
    std::ifstream f(file.c_str());
    
    if (!f.is_open())
    {
        std::cerr << "Could not open file " << file << " " << LINEINFO << std::endl;;
        exit(1);
        return;
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
        //x*= FMGEV;
        //y *= FMGEV;
        // Datafile step is 0.02 fm
        // Once we have load all points, we will sift all coordinates such that 0 is at the center
        
        double step = 0.04; // standard 128x128
	//double step = 0.0025;
	//double step=0.08;
        //double step = 0.007;
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
    

    
    // Now, given that we have a point (x,y), we can find the index xind such that
    // xcoords[xind] is closest to x, and similarly ind
    // Then, the corresponding Wilson line is
    // wilsonlines[ xcoords.size()*xind + yind]
    // Of course this is symmetric and we could just as well swap xind and yind

    
    std::cout <<"# Loaded " << wilsonlines.size() << " Wilson lines from file " << file << ", grid size " << xcoords.size() << " x " << ycoords.size() << " grid range [" << xcoords[0] << ", " << xcoords[xcoords.size()-1] << "]" << " step size " << xcoords[1]-xcoords[0] << " GeV^-1" << std::endl;

        
    
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


std::string IPGlasma::InfoStr()
{
    std::stringstream ss;
    ss << "# IPGlasma loaded from file " << datafile ;
    return ss.str();
}

std::vector<double> &IPGlasma::GetXCoordinates()
{
    return xcoords;
}
