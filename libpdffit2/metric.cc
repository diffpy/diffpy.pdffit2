#include <iostream>
#include <math.h>
#include "pdffit.h"

void dput(double a, double da)
{
	/*if (da != 0.0)
		cout << a << "(" << da << ") ";
	else
		cout << a << " ";*/
	
	cout << cc(a,da);
	
}

/*********************************************************************
    Calculates lattice constants, metric and reciprocal metric
     	tensor, permutation tensors and unit cell volume.
	It's done quite some old fashioned way, rather than calculating
	the direct metric tensor and its inverse.

	This extended version calculates standard deviations as well
**********************************************************************/
void Phase::lattice(bool lout)
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

	abscosa =abs(cosa);
	abscosb =abs(cosb);
	abscosg =abs(cosg);

	dcosa = abs(sina*rad*dwin[0]);
	dcosb = abs(sinb*rad*dwin[1]);
	dcosg = abs(sing*rad*dwin[2]);

	voll = 1.0 - sqr(cosa) - sqr(cosb) - sqr(cosg) + 2.0*cosa*cosb*cosg;

	if (voll > 0.0)
	{
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
			i1 = mod(i+1,3);    // i1 = 1,2,0
			i2 = mod(i+2,3);    // i2 = 2,0,1
			
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

//------ - output ?

		if (lout)
		{
			cout << " Phase number           : " << iphase << endl;
			cout << " Phase title            : " << name << endl << endl;

			cout << " Lattice parameters     : ";
			dput(a0[0] , da0[0]); dput(a0[1] , da0[1]); dput(a0[2] , da0[2]);
			cout << endl;
			
			cout << "           & angles     : ";
			dput(win[0] , dwin[0]); dput(win[1] , dwin[1]); dput(win[2] , dwin[2]);
			cout << endl;
				
			cout << " Unit cell volume       : ";
			dput(v, dv);
			cout << endl;

			cout << " Number density         : ";
			dput(rho0, drho0);
			cout << endl;
			
			for (j=0; j<3; j++)
			{
				if (!j) cout << " Metric tensor          : ";
				else    cout << "                          ";
				for (i=0; i<3; i++)
				{
					dput(gten[i][j], dgten[i][j]);
				}
				cout << endl;
			}
			cout << endl;

			cout << " Recip. lat. parameters : ";
			dput(ar[0] , dar[0]); dput(ar[1] , dar[1]); dput(ar[2] , dar[2]);
			cout << endl;

			cout << "               & angles : ";
			dput(wrez[0] , dwrez[0]); dput(wrez[1] , dwrez[1]); dput(wrez[2] , dwrez[2]);
			cout << endl;

			cout << " Recip. unit cell vol.  : ";
			dput(vr, dvr);
			cout << endl;

			for (j=0; j<3; j++)
			{
				if (!j) cout << " Recip. metric tensor   : ";
				else    cout << "                          ";
				for (i=0; i<3; i++)
				{
					dput(rten[i][j], drten[i][j]);
				}
				cout << endl;
			}
			cout << endl;
		}

	  // ????????  cr_la = .false.
	}
	else  // zero volume
		throw_exception(-35, ER_APPL);
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
				if (abs(ten[i][j]) < 1e-6) ten[i][j] = 0.0;
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
				if (dten[i][j] < 1e-6) dten[i][j] = 0.0;
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
double Phase::skalpro(double h[3], double k[3])
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
			dskalpro += abs(dh[i]*k[j]*gten[i][j]) + abs(h[i]*dk[j]*gten[i][j])
     	                  + abs(h[i]*k[j]*dgten[i][j]);
						  
	return dskalpro;
}

#if FORTRAN
c*****7*****************************************************************
	real function do_blen(lspace,u,v)
c-
c	Calculates the length of the vector v-u in the space defined by
c	lspace (.true. = real space , .false. = reciprocal space)
c+
	implicit       none
c
	include       'config.inc'
	include       'crystal.inc'
c
	logical        lspace
	real           u(3),v(3)
c
	integer        i
	real           w(3)
	real           skalpro
c
	do i=1,3
	  w(i) = v(i) - u(i)
	enddo
c
	if(lspace) then
	  do_blen = sqrt( skalpro(w,w,cr_gten,cr_ipha))
	else
	  do_blen = sqrt( skalpro(w,w,cr_rten,cr_ipha))
	endif
c
	end
c*****7*****************************************************************
	real function do_bang(lspace,u,v,w)
c-
c	Calculates the angle between vectors v-u and v-w in the space 
c	defined by lspace (.true. = real space , .false. = reciprocal space)
c+
	implicit       none
c
	include       'config.inc'
	include       'crystal.inc'
	include       'errlist.inc'
	include       'wink.inc'
c
	logical        lspace
	real           u(3),v(3),w(3)
c
	integer        i
	real           xx,xy,yy,arg
	real           x(3),y(3)
	real           skalpro
c
	do i=1,3
	  x(i) = v(i) - u(i)
	  y(i) = v(i) - w(i)
	enddo
c
	if(lspace) then
	  xx = sqrt( skalpro(x,x,cr_gten,cr_ipha))
	  xy =       skalpro(x,y,cr_gten,cr_ipha)
	  yy = sqrt( skalpro(y,y,cr_gten,cr_ipha))
	else
	  xx = sqrt( skalpro(x,x,cr_rten,cr_ipha))
	  xy =       skalpro(x,y,cr_rten,cr_ipha)
	  yy = sqrt( skalpro(y,y,cr_rten,cr_ipha))
	endif
	if(xx.gt.0 .and. yy.gt.0) then
	  arg = xy/(xx*yy)
	  if (abs(arg).le.1.0) then
	    do_bang = acos(arg)/rad
	  else
	    do_bang = 0.0
	  endif
	else
	  ier_num = -32
	  ier_typ = ER_APPL
	  do_bang = 0.0
	endif
c
	end
c*****7*********************************************************
c*****7*****************************************************************
c*****7*****************************************************************
	subroutine vekprod(u,v,ww,eps,rten,ipha)
c-
c	calculates the VECTORPRODUCT in general triclinic space
c	with  EPS and RTEN in direct space
c	with REPS and GTEN in reciprocal space
c+
	implicit       none
c
	include        'config.inc'
c
	integer        i,j,k,l,ipha
	real           u(3),v(3),ww(3)
	real           eps(3,3,3,MAXPHA),rten(3,3,MAXPHA)
c
	do i=1,3
	  ww(i)=0.0
	  do j=1,3
	    do k=1,3
	      do l=1,3
	        ww(i)=ww(i) + eps(j,k,l,ipha)*u(k)*v(l)*rten(j,i,ipha)
	      enddo
	    enddo
	  enddo
	enddo
c
	end
C*****7*****************************************************************
	subroutine matmulx(a,b,c)
c+
c	Matrixmultiplication c=a*b
c-
	implicit      none
c
	integer       i,j,k
	real          a(3,3),b(3,3),c(3,3)
c
	do i=1,3
	  do j=1,3
	    a(i,j)=0.0
	    do k=1,3
	      a(i,j)=a(i,j) + b(i,k)*c(k,j)
	    enddo
	  enddo
	enddo
	end
c*****7*****************************************************************
	subroutine transmat(mat,idim)
c-
c	Replaces a matrix by its transpose
c+
	implicit      none
c
	integer       idim,i,j
	real          mat(idim,idim),d
c
	do i=2,idim
	  do j=1,i-1
	    d        = mat(i,j)
	    mat(i,j) = mat(j,i)
	    mat(j,i) = d
	  enddo
	enddo
c
	end
c*****7*****************************************************************
	subroutine invmat(imat,a)
c
c	calculates the inverse matrix "imat" to input matrix "a"
c 
	implicit      none
c
	include      'errlist.inc'
c
	real imat(3,3),a(3,3),det
c
	det= a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) +
     &       a(2,1)*(a(3,2)*a(1,3) - a(1,2)*a(3,3)) +
     &       a(3,1)*(a(1,2)*a(2,3) - a(1,3)*a(2,2))

	if(abs(det).gt.0.0) then
	  ier_num = 0
	  ier_typ = ER_NONE
c
	  imat(1,1)= (a(2,2)*a(3,3) - a(2,3)*a(3,2) )/det
	  imat(1,2)=-(a(1,2)*a(3,3) - a(1,3)*a(3,2) )/det
	  imat(1,3)= (a(1,2)*a(2,3) - a(1,3)*a(2,2) )/det
c
	  imat(2,1)=-(a(2,1)*a(3,3) - a(3,1)*a(2,3) )/det
	  imat(2,2)= (a(1,1)*a(3,3) - a(1,3)*a(3,1) )/det
	  imat(2,3)=-(a(1,1)*a(2,3) - a(1,3)*a(2,1) )/det
c
	  imat(3,1)= (a(2,1)*a(3,2) - a(3,1)*a(2,2) )/det
	  imat(3,2)=-(a(1,1)*a(3,2) - a(1,2)*a(3,1) )/det
	  imat(3,3)= (a(1,1)*a(2,2) - a(1,2)*a(2,1) )/det
c
	else
	  ier_num = -1
	  ier_typ = ER_MATH
	endif
c
	end      
C*****7*****************************************************************
	subroutine invmat4(matrix)
c-
c	inverts a 4*4 Symmetry operation
c+
	implicit       none
c
	include        'errlist.inc'
c
	real           matrix(4,4)
	integer        i,j
	real           a(3,3),b(3,3),t(3)
c
	do i=1,3
	  do j=1,3
	    a(i,j) = matrix(i,j)
	  enddo
	  t(i) = matrix(i,4)
	enddo
c
	call invmat(b,a)
c
	if(ier_num.eq.0) then
	  do i=1,3
	    matrix(i,4) = 0.0
	    do j=1,3
	      matrix(i,j) = b(i,j)
	      matrix(i,4) = matrix(i,4) - b(i,j)*t(j)
	    enddo
	  enddo
	endif
c
	end
c*****7*****************************************************************
C*****7*****************************************************************
	subroutine trans(uc,gmat,ipha,up,idim)
c+
c	Transforms a point in the crystal space into plot space 
c	and vice versa
c-
	implicit	none
c
	include		'config.inc'
c
	integer	 i,j,idim,ipha
	real gmat(idim,idim,MAXPHA),uc(idim),up(idim)
c
	do i=1,idim
	  up(i)=0.0
	  do j=1,idim
	    up(i)=up(i)+gmat(i,j,ipha)*uc(j)
	  enddo
	enddo
c
	end
C*****7*****************************************************************
#endif
