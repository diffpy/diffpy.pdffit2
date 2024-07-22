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
* Custom output and error streams.
*
***********************************************************************/

#ifndef OUTPUTSTREAMS_H_INCLUDED
#define OUTPUTSTREAMS_H_INCLUDED

#include <iostream>

namespace NS_PDFFIT2 {

extern std::ostream* pout;
extern std::ostream* perr;

}

#endif	// OUTPUTSTREAMS_H_INCLUDED
