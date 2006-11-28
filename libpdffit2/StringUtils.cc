/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Utilities for string formatting.
*
* $Id$
*
***********************************************************************/

#include <cmath>
#include <sstream>
#include <iomanip>

#include "StringUtils.h"

using namespace std;

// return an uppercase copy of string s
string toupper(string s)
{
    for (string::iterator ii = s.begin(); ii != s.end(); ++ii)
    {
	*ii = toupper(*ii);
    }
    return s;
}

// Create format x(dx)
string putxdx(double x, double dx)
{
    ostringstream ostream;

    if (dx <= (1e-8*abs(x)) )
    {
	ostream << x;
    }
    else if (isnan(dx))
    {
	ostream << x << "(NaN)";
    }
    else 
    {
	const double rf=1-log10(9.5);  // rounding factor

	// compute exponents <ipowdx> and <ipowx> of dx and x
	// add (1-log10(9.5)) for <ipowdx> to allow for rounding of error

	int ipowdx = int(floor(log10(dx)+(dx<x ? rf : 0)));
	int ipowx =  int(floor(log10(abs(x))));

	// compute the exponent <base> of the standard deviation to transform
	// standard deviation in integer sd, i.e. nint(dx.mantissa)
	double base = exp10(ipowdx);
	int sd = int(rint(dx/base));  // 0.0025 -> 3


	// compute mantissa of x
	double mantissa = x/exp10(ipowx);

	// allow for rounding up of mantissa if dx has larger exponent than x 
	if ( (ipowdx > ipowx) && (mantissa >= 5) ) { ipowx++; mantissa = 1; }

	// Notation: 3 possibilities according to the value of exponent of dx and x: 
	//		edx <= 0, edx > 0 but edx <= ex and edx > 0 and edx > ex
	// (also added is a sophisticated test when edx <=0, using  
	// a lower exponent limit on dx and an upper exponent limit on x
	// to decide on using the scientific notation).
	if ( (ipowdx<=0) && (ipowdx>=-6) && (ipowx<=6) )
	    ostream << fixed << setprecision(-ipowdx) << x << "(" << sd << ")";
	else
	{
	    if (ipowx >= ipowdx)
	    {
		ostream << fixed << setprecision(ipowx-ipowdx) << mantissa << "(" << sd << ")";
		if (ipowx) ostream << "E" << showpos << ipowx << noshowpos;
	    }
	    else
	    {
		ostream << fixed << 0 << "(" << sd << ")";
		if (ipowdx) ostream << "E" << showpos << ipowdx << noshowpos;
	    }
	}
	ostream.unsetf(ios_base::fixed);  
	ostream << setprecision(6);
    }

    return ostream.str();
}

string cc(double x, double dx)
{
    ostringstream ostream;
	
    ostream << setw(15) << putxdx(x,dx);

    return ostream.str();
}

// End of file
