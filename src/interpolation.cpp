/*
 * AmplitudeLib tools
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2015
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation.hpp"

// This is also defined in config.hpp, but if this class is used standalone
// it is more safe to not to include config.hpp but define this here.
#ifndef LINEINFO
    #define LINEINFO __FILE__ << ":" << __LINE__
#endif

typedef unsigned int uint;
using std::isinf;
using std::isnan;

/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator::Initialize()
{
    int status=0;
    out_of_range_errors = true;
    if (ready)
    {
        // Interpolator is already initialized: clear it (GSL part)
        // and re-initialize
        switch(method)
        {
            case INTERPOLATE_SPLINE:
                if (spline != NULL)
                {
                    gsl_spline_free(spline);
                    spline=NULL;
                }
                if (acc != NULL)
                {
                    gsl_interp_accel_free(acc);
                    acc=NULL;
                }
                break;
            case INTERPOLATE_BSPLINE:
#ifdef ENABLE_BSPLINE
                gsl_bspline_free(bw);
                gsl_bspline_deriv_free(derbw);
                gsl_vector_free(B);
                gsl_matrix_free(X);
                gsl_vector_free(c);
                
                gsl_matrix_free(cov);
                gsl_multifit_linear_free(mw);
#endif
                
#ifndef ENABLE_BSPLINE
                std:cerr << "Bspline interpolation is not enabled! Recompile with ENABLE_BSPLINE!" << std::endl;
#endif
                break;
        }
        ready=false;
        
    }
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc(gsl_interp_cspline, points);
            status = gsl_spline_init(spline, xdata, ydata, points);
            break;
        case INTERPOLATE_BSPLINE:
#ifdef ENABLE_BSPLINE
            gsl_vector *x = gsl_vector_alloc(points);
            gsl_vector *y = gsl_vector_alloc(points);
            gsl_vector *w = gsl_vector_alloc(points);

            for (int i=0; i< points; i++)
            {
                gsl_vector_set(x, i, xdata[i]);
                gsl_vector_set(y, i, ydata[i]);
                gsl_vector_set(w, i, 1.0);
            }
     
            /* allocate a cubic bspline workspace (k = 4) */
            bw = gsl_bspline_alloc(k, nbreak);
            derbw = gsl_bspline_deriv_alloc(k);
            B = gsl_vector_alloc(ncoeffs);
       
            X = gsl_matrix_alloc(points, ncoeffs);
            c = gsl_vector_alloc(ncoeffs);
       
            cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
            mw = gsl_multifit_linear_alloc(points, ncoeffs);
     
     
            // use uniform breakpoints
            gsl_bspline_knots_uniform(xdata[0], xdata[points-1], bw);
     
            /* construct the fit matrix X */
            for (int i = 0; i < points; ++i)
            {
               double xi = gsl_vector_get(x, i);
             
               /* compute B_j(xi) for all j */
               gsl_bspline_eval(xi, B, bw);
             
               /* fill in row i of X */
               for (int j = 0; j < ncoeffs; ++j)
               {
                  double Bj = gsl_vector_get(B, j);
                  gsl_matrix_set(X, i, j, Bj);
               }
            }
     
            /* do the fit */
            double chisq;
            gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

            gsl_vector_free(x);
            gsl_vector_free(y);
            gsl_vector_free(w);
#endif
            break;
    }
    ready=true;
    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        return -1;
    }
    return 0;   //ok, there is no error handling at the moment...
}


double Interpolator::Evaluate(double x)
{
    if (isnan(x) or isinf(x))
    {
        cerr << "Trying to evaluate interpolator with x=" << x << " at " << LINEINFO << endl;
        exit(1);
    }
    
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator::Initialize()?" << endl;
        return 0;
    }

    if (x<minx or x>maxx)
    {
		if (freeze)
		{
			if (x<minx) return freeze_underflow;
			else return freeze_overflow;
		}
		if (x < 0.9999*minx or x > 1.00001*maxx)	// if not true, no need to display error
        {
            if (out_of_range_errors)
                cerr << "x=" << x << " is not within limits [" << minx << ", " << maxx << "], forcing "
                    << "it in that interval! " << LINEINFO << endl;
        }
        if (x<minx) x=minx*1.00001;
        if (x>maxx) x=maxx*0.999999;
    }
    
    double res, yerr; int status;
    res=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_e(spline, x, acc, &res);
            if (status)
            {
                cerr << "Interpolation failed at " << LINEINFO << ", error " << gsl_strerror(status)
                 << " (" << status << "), x=" << x << ", minx=" << xdata[0]
                 << ", maxx=" << xdata[points-1] << ", result=" << res << endl;
                 exit(1);
            }
            break;
        case INTERPOLATE_BSPLINE:
#ifdef ENABLE_BSPLINE
            gsl_bspline_eval(x, B, bw);
            gsl_multifit_linear_est(B, c, cov, &res, &yerr);

            /*if (std::abs(yerr/res)>0.05 )
            {
                cerr << "Interpolation failed at " << LINEINFO << ": bspline result "
                << res << " pm " << yerr << " relerr " << std::abs(yerr/res) << endl;
            }*/
#endif
            break;
        default:
            cerr << "Interpolation method is invalid! " << LINEINFO << endl;
            exit(1);
    }

    if (isnan(res) or isinf(res))
    {
        cerr << "Interpolation at x=" << x << " gives " << res << endl;
		return 0;
        exit(1);
    }

    
    return res;   
}

double Interpolator::Derivative(double x)
{
    double res=0; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
#ifdef ENABLE_BSPLINE
            gsl_matrix* mat = gsl_matrix_alloc(nbreak+k-2, 2);
            gsl_bspline_deriv_eval(x, 1, mat, bw, derbw);
            for (int i=0; i<ncoeffs; i++)
            {
                res += gsl_vector_get(c, i)*gsl_matrix_get(mat, i, 1);
            }
            gsl_matrix_free(mat);
#endif
            return res;
    }
    if (status)
        cerr << "An error occurred while evaluating the derivative at x=" << x
        << " result " << res << " " << LINEINFO << endl;

    return res;
}

double Interpolator::Derivative2(double x)
{
    double res; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv2_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
            cerr << "2nd derivative is not implemented for BSPLINE interpolation!"
            << " " << LINEINFO << endl;
            break;
    }

    if (status)
    {
        cerr << "2nd derivative interpolation failed at x=" << x <<", result "
        << res << " " << LINEINFO << endl;
    }
    return res;

}

Interpolator::Interpolator(double *x, double *y, int p)
{
    points=p;
    xdata=x;
    ydata=y;
    minx=x[0];
    maxx=x[p-1];
    method = INTERPOLATE_SPLINE;
    allocated_data=false;
    ready=false;
    freeze=false;
    freeze_underflow = y[0];
    freeze_overflow = y[p-1];

    for (int i=0; i<p; i++)
    {
        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
            }
        }
    }

    Initialize();
}

Interpolator::Interpolator(std::vector<double> &x, std::vector<double> &y)
{
    points = x.size();
    xdata = new double[points];
    ydata = new double[points];
    allocated_data=true;

    for (uint i=0; i<x.size(); i++)
    {
        xdata[i]=x[i];
        ydata[i]=y[i];

        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
                exit(1);
            }
        }
    }
    minx=xdata[0]; maxx=xdata[x.size()-1];
    method = INTERPOLATE_SPLINE;
    ready=false;
    freeze=false;
    freeze_overflow = y[y.size()-1];
    freeze_underflow = y[0];

    Initialize();
}

void Interpolator::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
    if (m == INTERPOLATE_BSPLINE)
        cerr << "BSPLINE interpolation should be tested in more detail before serious usage..." << endl;
}

void Interpolator::Clear()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            if (spline != NULL)
            {
                gsl_spline_free(spline);
                spline=NULL;
            }
            if (acc != NULL)
            {
                gsl_interp_accel_free(acc);
                acc=NULL;
            }
            break;
        case INTERPOLATE_BSPLINE:
#ifdef ENABLE_BSPLINE
            gsl_bspline_free(bw);
            gsl_bspline_deriv_free(derbw);
            gsl_vector_free(B);
            gsl_matrix_free(X);
            gsl_vector_free(c);
            
            gsl_matrix_free(cov);
            gsl_multifit_linear_free(mw);
#endif
            break;
    }

    if (allocated_data)
    {
        delete[] xdata;
        delete[] ydata;
        allocated_data=false;
    }
}

Interpolator::~Interpolator()
{
    Clear();

}

double* Interpolator::GetXData()
{
    return xdata;
}
double* Interpolator::GetYData()
{
    return ydata;
}
int Interpolator::GetNumOfPoints() const
{
    return points;
}
INTERPOLATION_METHOD Interpolator::GetMethod() const
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator::Interpolator(const Interpolator& inter)
{
    points=inter.GetNumOfPoints();
    xdata = new double[points];
    ydata = new double[points];
    allocated_data=true;


    gsl_spline *tmpspline = inter.GetGslSpline();
    for (int i=0; i<points; i++)
    {
        xdata[i] = tmpspline->x[i];
        ydata[i] = tmpspline->y[i];
    }
    minx = xdata[0]; maxx=xdata[points-1];
    method = inter.GetMethod();
    ready=false;
    Initialize();
}

gsl_spline* Interpolator::GetGslSpline() const
{
    return spline;
}

double Interpolator::MinX()
{
	return minx;
}

double Interpolator::MaxX()
{
	return maxx;
}


bool Interpolator::Freeze()
{
	return freeze;
}
void Interpolator::SetFreeze(bool f)
{
	freeze=f;
}
void Interpolator::SetUnderflow(double min)
{
	freeze_underflow=min;
}
 void Interpolator::SetOverflow(double max)
 {
	 freeze_overflow=max;
 }
double Interpolator::UnderFlow()
{
	 return freeze_underflow;
}
double Interpolator::OverFlow()
{
	return freeze_overflow;
}

void Interpolator::SetMaxX(double x)
{
    maxx=x;
}

void Interpolator::SetMinX(double x)
{
    minx=x;
}

void Interpolator::SetOutOfRangeErrors(bool er)
{
    out_of_range_errors=er;
}
