/*
 * General purpose interpolation class
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H



#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;


/**
 * Interpolation method, spline goes trough every
 * datapoint, bspline fits a noisy dataset.
 */
enum INTERPOLATION_METHOD {
    INTERPOLATE_SPLINE,
    INTERPOLATE_BSPLINE
};

// Enable support for bspline interpolation, requires gsl 1.x, not gsl 2
// #define ENABLE_BSPLINE

/**
 * Interpolates given data using spline (goes trough every data point)
 * or bspline (=noisy data) using the GSL routines.
 */
class Interpolator
{
    public:
        /**
         * Create interpolator from two arrays.
         *
         * Pointers to the arrays are saved and used when Interpolator::Initialize()
         * is called. Also these pointers can be asked from the class later (see GetXData()
         * and GetYData()). The arrays are not used for evaluating the interpolator.
         * @param x array of x coordinates
         * @param y array of y coordinates
         * @param p number of points in arrays
         */
        Interpolator(double* x, double* y, int p);

        /**
         * Create interpolator from two std::vectors.
         *
         * The given vectors are not saved or referenced later.
         */
        Interpolator(std::vector<double> &x, std::vector<double> &y);
        Interpolator(const Interpolator& inter);
        ~Interpolator();
        void Clear();
        /**
         * Evaluate interpolator f(x)
         */
        double Evaluate(double x);
        /**
         * Evaluate 1st derivative of the interpolated function f'(x)
         */
        double Derivative(double x);    
        /**
         * Evaluate 2nd derivative of the interpolated function f''(x)
         */
        double Derivative2(double x);   
        /**
         * Select interpolation method (spline or bspline)
         */
        void SetMethod(INTERPOLATION_METHOD m);

        /**
         * Initialize interpolator.
         *
         * This is done automatically at the constructor. If
         * interpolation method is changed (see SetMethod()), this
         * must be evaluated.
         */
        int Initialize();

        /**
         * Largest x value supported
         */
        double MinX();
        /**
         * Smallest x-value supported
         */
        double MaxX();

        double* GetXData();
        double* GetYData();
        gsl_spline* GetGslSpline() const;
        int GetNumOfPoints() const;
        INTERPOLATION_METHOD GetMethod() const;
        
        bool Freeze();

        /**
         * Set what to do when interpolated function is evaluated
         * outside the interpolation region.
         *
         * @param f true if pre-saved values are returned
         */
        void SetFreeze(bool f);
        /**
         * Set value that is returned if interpolator is evaluated at
         * x<minx
         */
        void SetUnderflow(double min);
        /**
         * Set value that is returned if interpolator is evaluated at
         * x>maxx
         */
        void SetOverflow(double max);
        double UnderFlow();
        double OverFlow();

        /**
         * Set interpolation range (maxx), default is upper limit of the
         * xdata
         */
        void SetMaxX(double x);
        /**
         * Set interpolation range (minx), default is lower  limit of the
         * xdata
         */
        void SetMinX(double x);

        /**
         * Choose wether an error will be printed to stderr if interpolator
         * is evaluated outside the interpolation range.
         */
        void SetOutOfRangeErrors(bool set);

    private:
        INTERPOLATION_METHOD method;
        double* xdata, *ydata;
        bool allocated_data;    // true if we allocated xdata,ydata
        double minx,maxx;
        int points;
        bool ready;
        
        bool freeze;		// true if return freeze_under/overflow if
        double freeze_underflow;	// asked to evaluate interpolation
		double freeze_overflow;	// outside the spesified range
        
        // spline
        gsl_interp_accel *acc;
        gsl_spline *spline;
#ifdef ENABLE_BSPLINE
        // bspline
        gsl_bspline_workspace *bw;
        gsl_bspline_deriv_workspace *derbw;
        gsl_vector *B;
        gsl_vector *c;
        gsl_matrix *X;
        gsl_matrix *cov;
        gsl_multifit_linear_workspace *mw;
#endif
        static const int k=4;
        static const int ncoeffs = 12;
        static const int nbreak = ncoeffs-k+2;

        bool out_of_range_errors;


};




#endif
