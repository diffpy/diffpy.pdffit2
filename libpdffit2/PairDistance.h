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
* Simple struct for pair distance and indices of atoms from the pair
*
* $Id$
*
***********************************************************************/

#ifndef PAIRDISTANCE_H_INCLUDED
#define PAIRDISTANCE_H_INCLUDED

struct PairDistance
{
    // Data members
    double dij;
    double ddij;
    int i;
    int j;

    // Constructor
    PairDistance() : dij(0.0), ddij(0.0), i(0), j(0)
    { }

};

// Comparison operator
inline bool operator<(const PairDistance& pd0, const PairDistance& pd1)
{
    return (pd0.dij < pd1.dij);
}

#endif	// PAIRDISTANCE_H_INCLUDED
