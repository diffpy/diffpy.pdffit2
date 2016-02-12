/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Utilities from numerical recipies.
*
* Comments:
*
***********************************************************************/

#ifndef NRUTIL_H_INCLUDED
#define NRUTIL_H_INCLUDED

#include <cmath>

void nrerror(char error_text[]);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

#endif	// NRUTIL_H_INCLUDED
