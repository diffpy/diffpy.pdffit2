/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Math functions and numerical constants for pdffit2 and formula parser.
*
* Comments:
*
***********************************************************************/

#ifndef MATHUTILS_H_INCLUDED
#define MATHUTILS_H_INCLUDED

#include <cmath>
#include <limits>

// numerical constants
const double rad = M_PI/180.0;
const double double_eps = (1.0 + sqrt(std::numeric_limits<double>().epsilon())) - 1.0;
const double double_max = std::numeric_limits<double>().max();
const double deltar_tol = 1.0e-3;

// MS compatibility fix - define functions not provided by MSVC cmath
#ifdef _MSC_VER
#include <cfloat>

#if _MSC_VER < 1800 // log2 has been available since MSVC 2013
inline double log2(double x)	{ return log(x)/log(2.0); }
#endif

inline int isnan(double x)	{ return _isnan(x); }
inline double round(double x)	{ return (x < 0) ? ceil(x - 0.5) : floor(x + 0.5); }

#endif	// _MSC_VER

// nearest integer
inline int nint(const double x)
{
    return (int) round(x);
}

// math functions and derivatives used in formula parser

inline double neg(double x)	{ return -x; }
inline double dneg(double x)	{ return -1.0; }

inline double sqr(double x)	{ return x*x; }
inline double dsqr(double x)	{ return 2.0*x; }

inline double cube(double x)	{ return x*x*x; }
inline double dcube(double x)	{ return 3.0*sqr(x); }

inline double dsqrt(double x)	{ return 0.5/sqrt(x); }

inline double dexp(double x)	{ return exp(x); }
inline double dlog(double x)	{ return 1.0/x; }

inline double sind(double x)	{ return sin(rad*x); }
inline double dsind(double x)	{ return rad*cos(rad*x); }

inline double cosd(double x)	{ return cos(rad*x); }
inline double dcosd(double x)	{ return -rad*sin(rad*x); }

inline double tand(double x)	{ return tan(rad*x); }
inline double dtand(double x)	{ return rad/sqr(cosd(x)); }

inline double dsin(double x)	{ return cos(x); }
inline double dcos(double x)	{ return -sin(x); }
inline double dtan(double x)	{ return 1.0/sqr(cos(x)); }

inline double dasin(double x)	{ return 1.0/sqrt(1.0 - x*x); }
inline double dacos(double x)	{ return -1.0/sqrt(1.0 - x*x); }
inline double datan(double x)	{ return 1/(1+sqr(x)); }

inline double asind(double x)	{ return asin(x)/rad; }
inline double dasind(double x)	{ return dasin(x)/rad; }

inline double acosd(double x)	{ return acos(x)/rad; }
inline double dacosd(double x)	{ return dacos(x)/rad; }

inline double atand(double x)	{ return atan(x)/rad; }
inline double datand(double x)	{ return datan(x)/rad; }

#endif	// MATHUTILS_H_INCLUDED
