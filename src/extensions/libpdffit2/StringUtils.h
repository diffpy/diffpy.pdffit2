/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Utilities for string formatting.
*
***********************************************************************/

#ifndef STRINGUTILS_H_INCLUDED
#define STRINGUTILS_H_INCLUDED

#include <string>

std::string toupper(std::string s);

class FormatValueWithStd
{
    private:

	// Data members
	int f_width;
	bool f_left;
	bool f_leading_blank;
	int f_std_precision;

    public:

	// Constructor:
	FormatValueWithStd()
	{
	    f_width = 0;
	    f_left = false;
	    f_leading_blank = false;
	    f_std_precision = 2;
	}

	// Methods:
	std::string operator() (double x, double dx);
	inline FormatValueWithStd& width(int w)
	{
	    f_width = w;
	    return *this;
	}
	inline FormatValueWithStd& left()
	{
	    f_left = true;
	    return *this;
	}
	inline FormatValueWithStd& right()
	{
	    f_left = false;
	    return *this;
	}
	inline FormatValueWithStd& leading_blank(bool flag)
	{
	    f_leading_blank = flag;
	    return *this;
	}
};

#endif	// STRINGUTILS_H_INCLUDED
