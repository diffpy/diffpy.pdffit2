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

#include "OutputStreams.h"

std::ostream* NS_PDFFIT2::pout = &std::cout;
std::ostream* NS_PDFFIT2::perr = &std::cerr;

// End of file
