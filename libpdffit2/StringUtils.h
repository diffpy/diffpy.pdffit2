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
* $Id$
*
***********************************************************************/

#ifndef STRINGUTILS_H_INCLUDED
#define STRINGUTILS_H_INCLUDED

#include <string>

std::string toupper(std::string s);
std::string putxdx(double x, double dx);
std::string cc(double x, double dx);

#endif	// STRINGUTILS_H_INCLUDED
