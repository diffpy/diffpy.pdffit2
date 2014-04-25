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
* Definitions of several Phase methods
*
* Comments:
*
***********************************************************************/

#include <iostream>

#include "MathUtils.h"
#include "StringUtils.h"
#include "pdffit.h"
using NS_PDFFIT2::pout;

/*********************************************************************
  Calculates lattice constants, metric and reciprocal metric
  tensor, permutation tensors and unit cell volume.
  It's done quite some old fashioned way, rather than calculating
  the direct metric tensor and its inverse.

  This extended version calculates standard deviations as well
 **********************************************************************/
void Phase::lattice()
{
    double abscosa, abscosb, abscosg, arg, darg;
    double dcosa, dcosb, dcosg, voll, dvoll;
    double cos1, cos2, cosi, sin1, sin2, sini;
    double dcos1, dcos2, dcosi, dsin1, dsin2, dsini;
    int i, i1, i2, j, k, ncc;

    cosa = cosd(win[0]);
    cosb = cosd(win[1]);
    cosg = cosd(win[2]);

    sina = sind(win[0]);
    sinb = sind(win[1]);
    sing = sind(win[2]);

    abscosa = fabs(cosa);
    abscosb = fabs(cosb);
    abscosg = fabs(cosg);

    dcosa = fabs(sina*rad*dwin[0]);
    dcosb = fabs(sinb*rad*dwin[1]);
    dcosg = fabs(sing*rad*dwin[2]);

    voll = 1.0 - sqr(cosa) - sqr(cosb) - sqr(cosg) + 2.0*cosa*cosb*cosg;

    if (voll <= double_eps)
    {
	throw structureError("Unit cell volume is not positive.");
    }
    v = sqrt(voll)*a0[0]*a0[1]*a0[2];

    dvoll = 2.0*(abscosa*dcosa + abscosb*dcosb + abscosg*dcosg +
	    dcosa*abscosb*abscosg + abscosa*dcosb*abscosg + abscosa*abscosb*dcosg);
    dv = 0.5/sqrt(voll)*dvoll*a0[0]*a0[1]*a0[2]
	+ sqrt(voll)*(da0[0]*a0[1]*a0[2] + a0[0]*da0[1]*a0[2] + a0[0]*a0[1]*da0[2]);

    vr = 1.0/v;
    dvr = 1.0/(v*v)*dv;

    //------ - calculate direct metric tensor

    tensor(gten,a0,win);
    dtensor(a0,win,dgten,da0,dwin);

    //------ - calculate reciprocal lattice constants

    // WHY ARE THERE NO ABSOLUTE VALUES IN THE FOLLOWING FORMULAE FOR THE STANDARD
    // DEVIATIONS??
    for (i=0; i<3; i++)     // i = 0,1,2
    {
	i1 = (i+1) % 3;    // i1 = 1,2,0
	i2 = (i+2) % 3;    // i2 = 2,0,1

	cos1 = cosd(win[i1]);
	cos2 = cosd(win[i2]);
	cosi = cosd(win[i]);
	sin1 = sind(win[i1]);
	sin2 = sind(win[i2]);
	sini = sind(win[i]);

	ar[i] = a0[i1]*a0[i2]*sini/v;
	arg = (cos1*cos2-cosi)/(sin1*sin2);
	wrez[i] = acosd(arg);

	dcos1 = sin1*rad*dwin[i1];
	dcos2 = sin2*rad*dwin[i2];
	dcosi = sini*rad*dwin[i];
	dsin1 = cos1*rad*dwin[i1];
	dsin2 = cos2*rad*dwin[i2];
	dsini = cosi*rad*dwin[i];

	dar[i] = ( (da0[i1]*a0[i2] + a0[i1]*da0[i2])*sini +
		a0[i1]*a0[i2]*dsini) /v + a0[i1]*a0[i2]*sini/(v*v)*dv;
	darg = (dcos1*cos2 + cos1*dcos2 + dcosi)/(sin1*sin2)
	    + arg/sin1*dsin1 + arg/sin2*dsin2;
	dwrez[i] = 1./sqrt(1-arg*arg)*darg/rad;
    }

    //------ - calculate reciprocal tensor

    tensor(rten,ar,wrez);
    dtensor(ar,wrez,drten,dar,dwrez);

    //------ - calculate permutation tensors

    for (i=0; i<3; i++)
    {
	for (j=0; j<3; j++)
	{
	    for (k=0; k<3; k++)
	    {
		eps(i,j,k)   = 0.0;
		reps(i,j,k)  = 0.0;
		deps(i,j,k)  = 0.0;
		dreps(i,j,k) = 0.0;
	    }
	}
    }

    eps(0,1,2)   =  v;
    eps(1,2,0)   =  v;
    eps(2,0,1)   =  v;
    eps(0,2,1)   = -v;
    eps(2,1,0)   = -v;
    eps(1,0,2)   = -v;
    reps(0,1,2)  =  vr;
    reps(1,2,0)  =  vr;
    reps(2,0,1)  =  vr;
    reps(0,2,1)  = -vr;
    reps(2,1,0)  = -vr;
    reps(1,0,2)  = -vr;

    deps(0,1,2)  =  dv;
    deps(1,2,0)  =  dv;
    deps(2,0,1)  =  dv;
    deps(0,2,1)  = -dv;
    deps(2,1,0)  = -dv;
    deps(1,0,2)  = -dv;
    dreps(0,1,2) =  dvr;
    dreps(1,2,0) =  dvr;
    dreps(2,0,1) =  dvr;
    dreps(0,2,1) = -dvr;
    dreps(2,1,0) = -dvr;
    dreps(1,0,2) = -dvr;

    //------ - Calculate number density

    np = 0.0;
    dnp = 0.0;

    for (i=0; i<natoms; i++)
    {
	np += atom[i].occ;
	dnp += atom[i].docc;
    }
    ncc = icc[0]*icc[1]*icc[2];
    rho0 = np/(ncc*v);

    //!!!! 	  error on rho0 seems very FISHY !!!
    drho0 = sqr(dnp/v) + sqr(dv*np/sqr(v));
    drho0 = sqrt(drho0);

}

void Phase::show_lattice()
{
    FormatValueWithStd value_std;
    value_std.leading_blank(true).left();

    *pout << " Phase number           : " << iphase << endl;
    *pout << " Phase title            : " << name << endl << endl;

    value_std.leading_blank(true).left();
    *pout << " Lattice parameters     :"
	<< value_std.width(20)(a0[0], da0[0])
	<< value_std.width(20)(a0[1], da0[1])
	<< value_std.width(0)(a0[2], da0[2]) << '\n';

    *pout << "           & angles     :"
	<< value_std.width(20)(win[0], dwin[0])
	<< value_std.width(20)(win[1], dwin[1])
	<< value_std.width(0)(win[2], dwin[2]) << '\n';

    *pout << " Unit cell volume       :"
	<< value_std.width(0)(v, dv) << endl;

    *pout << " Number density         :"
	<< value_std(rho0, drho0) << endl;

    for (size_t j = 0; j != 3; j++)
    {
	if (!j) *pout << " Metric tensor          :";
	else    *pout << "                         ";
	for (size_t i = 0; i != 3; i++)
	{
	    *pout << value_std.width(i < 2 ? 20 : 0)(gten[i][j], dgten[i][j]);
	}
	*pout << endl;
    }
    *pout << endl;

    *pout << " Recip. lat. parameters :"
	<< value_std.width(20)(ar[0], dar[0])
	<< value_std.width(20)(ar[1], dar[1])
	<< value_std.width(0)(ar[2], dar[2]) << endl;

    *pout << "               & angles :"
	<< value_std.width(20)(wrez[0], dwrez[0])
	<< value_std.width(20)(wrez[1], dwrez[1])
	<< value_std.width(0)(wrez[2], dwrez[2]) << endl;

    *pout << " Recip. unit cell vol.  :"
	<< value_std(vr, dvr) << endl;

    for (size_t j = 0; j != 3; j++)
    {
	if (!j) *pout << " Recip. metric tensor   :";
	else    *pout << "                         ";
	for (size_t i = 0; i != 3; i++)
	{
	    *pout << value_std.width(i < 2 ? 20 : 0)(rten[i][j], drten[i][j]);
	}
	*pout << '\n';
    }
    *pout << endl;
}



/*********************************************************************
  Calculates the metric tensor. Works both for direct and
  reciprocal metric tensor.
 **********************************************************************/
void Phase::tensor(double ten[3][3], double vec[3], double win[3])
{
    //include		'config.inc'

    const int idim=3;

    int i, j;

    for (i=0; i<idim; i++)
    {
	for (j=0; j<idim; j++)
	{
	    if(i != j)
	    {
		ten[i][j] = vec[i]*vec[j]*cosd(win[3-(i+j)]);
		if (fabs(ten[i][j]) < double_eps) ten[i][j] = 0.0;
	    }
	    else
		ten[i][j] = vec[i]*vec[j];
	}
    }
}

/***********************************************************************
  Calculates the standard devations for the metric tensor.
  Works both for direct and reciprocal metric tensor.
 ***********************************************************************/
void Phase::dtensor(double vec[3], double win[3], double dten[3][3],
	double dvec[3], double dwin[3])
{
    //include		'config.inc'
    //include		'wink.inc'

    const int idim=3;

    int	i, j;

    for (i=0; i<idim; i++)
    {
	for (j=0; j<idim; j++)
	{
	    if(i != j)
	    {
		dten[i][j] = dvec[i]*vec[j]*cosd(win[3-(i+j)]) +
		    vec[i]*dvec[j]*cosd(win[3-(i+j)])+
		    vec[i]*vec[j]*sind(win[3-(i+j)])*rad*dwin[3-(i+j)];
		if (dten[i][j] < double_eps) dten[i][j] = 0.0;
	    }
	    else
		dten[i][j] = dvec[i]*vec[j] + vec[i]*dvec[j];
	}
    }
}

/********************************************************
  Calulates the SCALARPRODUCT of two vectors
  1/D**2 = H(I)*K(J)*RTEN(I,J)
  uses the phase's metric
 *********************************************************/
double Phase::skalpro(const double h[3], const double k[3])
{
    const int idim=3;
    int i,j;
    double skalpro;

    skalpro = 0.0;
    for (i=0; i<idim; i++)
	for (j=0; j<idim; j++)
	    skalpro += h[i] * k[j] * gten[i][j];

    return skalpro;
}

/******************************************
  c	Calculates the error of scalar product
 *******************************************/
double Phase::dskalpro(double h[3] ,double k[3], double dh[3], double dk[3])
{
    const int idim=3;
    int	i, j;
    double dskalpro;

    dskalpro=0.0;
    for (i=0; i<idim; i++)
	for (j=0; j<idim; j++)
	    dskalpro += fabs(dh[i]*k[j]*gten[i][j]) + fabs(h[i]*dk[j]*gten[i][j])
		+ fabs(h[i]*k[j]*dgten[i][j]);

    return dskalpro;
}

// End of file
