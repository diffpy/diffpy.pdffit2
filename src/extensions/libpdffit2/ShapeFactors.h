/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* PDF corrections due to particle shape.
*
* Comments:
*    PDF correction for spherical particles obtained from
*    R.C. Howell et al, Phys. Rev. B 73, 094107 (2006)
*    http://link.aps.org/abstract/PRB/v73/e094107
*
***********************************************************************/

#ifndef SHAPEFACTORS_H_INCLUDED
#define SHAPEFACTORS_H_INCLUDED

#include <cmath>

// envelope function for spherical nanoparticle with diameter d
inline double sphereEnvelope(double r, double d)
{
    double rdratio = r/d;
    return (rdratio < 1.0) ? (1.0 - 1.5*rdratio + 0.5*pow(rdratio, 3)) : 0.0;
}

// derivative of sphereEnvelope by diameter d
inline double dsphereEnvelope(double r, double d)
{
    double rdratio = r/d;
    return (rdratio < 1.0) ? (1.5*rdratio/d - 1.5*pow(rdratio, 3)/d) : 0.0;
}

#endif	// SHAPEFACTORS_H_INCLUDED
