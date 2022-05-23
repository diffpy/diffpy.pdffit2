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
***********************************************************************/

#include <sstream>
#include <iomanip>

#include "StringUtils.h"
#include "MathUtils.h"

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

// class FormatValueWithStd
string FormatValueWithStd::operator() (double x, double dx)
{
    ostringstream out(f_leading_blank ? " " : "", ios_base::app);
    out << x << setprecision(f_std_precision);
    // do not write dx when it is too small
    if (dx > fabs(x)*1e-8)	out << " (" << dx << ')';
    else if (std::isnan(dx))		out << " (NaN)";
    // left-pad string to the width
    string rv = out.str();
    int rvlen = rv.size();
    // pad or prepend blanks as necessary
    if (rvlen < f_width)
    {
        size_t nblanks = f_width - rvlen;
        if (f_left) rv.append(nblanks, ' ');
        else	    rv.insert(0, nblanks, ' ');
    }
    return rv;
}

// End of file
