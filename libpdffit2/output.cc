/***********************************************************************
c	This file contains routines to display information.
c
c	Version: 1.0
c	Date:    01 June 1999
c	Author:  Th. Proffen (tproffen@lanl.gov)
**********************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include "pdffit.h"
using namespace std;
#include <iomanip>



/*********************************************
c	This routine creates the format x(dx)
*********************************************/
string putxdx(double x, double dx)
{
	ostringstream ostream;
	
	if (dx <= (1e-8*abs(x)) )
	{
		ostream << x;
	}
	else 
	{
		const double rf=1-log10(9.5);  // rounding factor
		
		// compute exponents <ipowdx> and <ipowx> of dx and x
		// add (1-log10(9.5)) for <ipowdx> to allow for rounding of error
		
		int ipowdx = int(floor(log10(dx)+(dx<x ? rf : 0)));
		int ipowx =  int(floor(log10(abs(x))));

		// compute the exponent <base> of the standard deviation to transform
		// standard deviation in integer sd, i.e. nint(dx.mantissa)
		double base = exp10(ipowdx);
		int sd = nint(dx/base);  // 0.0025 -> 3
		
		
		// compute mantissa of x
		double mantissa = x/exp10(ipowx);
		
		// allow for rounding up of mantissa if dx has larger exponent than x 
		if ( (ipowdx > ipowx) && (mantissa >= 5) ) { ipowx++; mantissa = 1; }

		//cout << ipowx << " " << ipowdx << " ";
		
		// Notation: 3 possibilities according to the value of exponent of dx and x: 
		//		edx <= 0, edx > 0 but edx <= ex and edx > 0 and edx > ex
		// (also added is a sophisticated test when edx <=0, using  
		// a lower exponent limit on dx and an upper exponent limit on x
		// to decide on using the scientific notation).
		if ( (ipowdx<=0) && (ipowdx>=-6) && (ipowx<=6) )
			ostream << fixed << setprecision(-ipowdx) << x << "(" << sd << ")";
		else
		{
			if (ipowx >= ipowdx)
			{
				ostream << fixed << setprecision(ipowx-ipowdx) << mantissa << "(" << sd << ")";
				if (ipowx) ostream << "E" << showpos << ipowx << noshowpos;
			}
			else
			{
				ostream << fixed << 0 << "(" << sd << ")";
				if (ipowdx) ostream << "E" << showpos << ipowdx << noshowpos;
			}
		}
		ostream.unsetf(ios_base::fixed);  
		ostream << setprecision(6);
	}

	return ostream.str();
}

string cc(double x, double dx)
{
	ostringstream ostream;
	
	ostream << setw(15) << putxdx(x,dx);

	return ostream.str();
}

static void test()
{
	/*_pp(cout.precision());
	cout << 1234.56789 << " " << 123.123456 << endl;
	cout << scientific << 1234.56789 << endl;
	cout << fixed << 1234.56789 << endl;
	cout.unsetf( ios_base::fixed );  
	cout << 1.0 << " " << 123.123456 << endl;*/
	
	//  Some # to test the use of cc(x,dx)
	cout << cc(1.2,0.0025) << endl;
	cout << cc(1.2,0.00256) << endl;
	cout << cc(1.2126784,0.00256) << endl;
	cout << cc(162.0415,0.00256) << endl;
	cout << cc(0.041,0.00256) << endl;
	cout << cc(0.0041,0.00256) << endl;
	cout << cc(0.00051,0.00256) << endl;
	cout << cc(0.00041,0.00256) << endl;
	cout << cc(234895.6, 0.001) << endl;
	cout << cc(95.6432, 0.096) << endl;	
	cout << cc(234895.612, 0.1) << endl;
	cout << cc(1.2346e-16, 1.456e-24) << endl;
	cout << cc(1.2346e12, 0.012) << endl;
	cout << cc(1.2346e-16, 0.012) << endl;
	cout << endl;
	
	cout << cc(351.6, 2.465) << endl;
	cout << cc(351, 0.0004) << endl;
	cout << cc(234895.6, 43.465) << endl;
	cout << cc(3456.34244,7.66) << endl;
	cout << cc(3456.34244,9.56) << endl;
	cout << cc(234.356,1.234) << endl;
	cout << cc(252.6, 13.465) << endl;
	cout << cc(252.6, 134.65) << endl;
	cout << cc(2.6, 13.465) << endl;
	cout << cc(9.6, 13.465) << endl;
	cout << cc(0.026, 13.465) << endl;
	cout << cc(1000000000000.0, 456100) << endl;
	cout << cc(1e12, 4.561e5) << endl;

	cout << cc(-1223.455,12.3) << endl;
}

/************************************************************
c	Outputs phase information for phase 'ipha' on file 'id'
    Thu Oct 13 2005 - CLF
    Modified code to handle general types of
    streams.
*************************************************************/
void Phase::output(ostream &fout)
{
	int ia, i;
	
	//test();
	fout << " " << string(78,'-') << endl
		<< " PHASE " << iphase << " : " << name << endl
		<< " " << string(78,'-') << endl;

	fout << " Scale factor          : " << cc(skal,dskal) << endl;

	if (corr_max > 0.0)
		fout << " Correlation limit [A] : " << corr_max << endl;
	
	fout << " Quad. corr. factor    : " << cc(delta,ddelta) << endl;

	fout << " Lin. corr. factor     : " << cc(gamma,dgamma) << endl;

	fout << " Low r sigma ratio     : " << cc(srat,dsrat) << endl
    	 << " R cutoff [A]          : " << cc(rcut,0.0) << endl;

	fout << " Lattice parameters    : ";
	for (i=0; i<3; i++)
		fout << cc(a0[i], da0[i]);
	fout << endl;

	fout << "           & angles    : ";
	for (i=0; i<3; i++)
		fout << cc(win[i], dwin[i]);
	fout << endl << endl;

	fout << " Atom positions & occupancies : " << endl;
	for (ia=0; ia<natoms; ia++)
	{
		Atom &atom=this->atom[ia];
		
		fout << "   " << setw(4) << at_lis[atom.iscat];
		
		for (i=0; i<3; i++)
			fout << cc(atom.pos[i],atom.dpos[i]);

		fout << cc(atom.occ,atom.docc) << endl;
	}
	fout << endl;

	fout << " Anisotropic temperature factors : " << endl;

	for (ia=0; ia<natoms; ia++)
	{
		Atom &atom=this->atom[ia];

		fout << "   " << setw(4) << at_lis[atom.iscat];
		
		for (i=0; i<3; i++)
			fout << cc(atom.u[i],atom.du[i]);
		fout << endl;

		if ( atom.u[3] || atom.u[4] || atom.u[5])
		{
			for (i=3; i<6; i++)
				fout << "            " << cc(atom.u[i],atom.du[i]);
			fout << endl;
		}
	}
}

/***********************************************************
c	Outputs information about this data set 
    Thu Oct 13 2005 - CLF
    Modified code to handle general types of
    streams.
************************************************************/
void DataSet::output(ostream &fout)
{
	fout << " " << string(78,'-') << endl
		 << " DATA SET : " << iset << " (" << name << ")" << endl
		 << " " << string(78,'-') << endl;

	fout << " Data range in r [A]   : " << setw(8) << rmin << " -> " << setw(8) << rmax
     	 << "      Step dr  : " << setw(8) << deltar << endl;
	
	fout << " Calculated range      : " << setw(8) << rcmin << " -> " << setw(8) << rcmax << endl;

	fout << " Refinement r range    : " << setw(8) << rfmin << " -> " << setw(8) << rfmax 
		 << "      Data pts : " << setw(5) << nfmin << " -> " << setw(5) << nfmax << endl;

	//if (lref) 
	//	fout << " Reference PDF file    : " << rname;

	fout << endl << " Experimental settings : " << endl;

	if (lxray)
		fout << "   Radiation           : X-Rays\n";
	else
		fout << "   Radiation           : Neutrons\n";

	if (qmax == 0.0)
		fout << "   Termination at Qmax : not applied\n";
	else
		fout << "   Termination at Qmax : " << setw(15) << qmax << " A**-1\n";

	if (sigmaq == 0.0)
		fout << "   DQ dampening Qsig   : not applied\n";
	else
		fout << "   DQ dampening Qsig   : " << cc(sigmaq,dsigmaq) << " A**-1\n";

	if (qalp == 0.0)
		fout << "   DQ broadening Qalp  : not applied\n";
	else
		fout << "   DQ broadening Qalp  : " << cc(qalp,dqalp) << " A**-1\n";

	fout << "   Scale factor        : " << cc(skal,dskal) << endl;

	fout << endl;

	fout << " Selected phases and atoms for this data set :\n";

	for(unsigned int ip=0; ip<psel.size(); ip++)
	{
		
		if (psel[ip])
		{
			fout << "   Phase " << ip+1 << " :\n";
		
			fout << "     Atoms (i) :  ";
			for(int i=0; i<psel[ip]->nscat; i++)
			{
				if (allowed_i[ip][i])
					fout << psel[ip]->at_lis[i] << "  ";
			}
			fout << endl;
			
			fout << "     Atoms (j) :  ";
			for(int i=0; i<psel[ip]->nscat; i++)
			{
				if (allowed_j[ip][i])
					fout << psel[ip]->at_lis[i] << "  ";
			}
			fout << endl;
		}
	}
}

/*************************************************
c	Outputs parameter information
    Thu Oct 13 2005 - CLF
    Modified code to handle general types of
    streams.
**************************************************/
void Fit::output(ostream &fout)
{

	fout << " " << string(78,'-') << endl
		 << " PARAMETER INFORMATION :" << endl
		 << " " << string(78,'-') << endl;

	int npar = 0;

	for (int i=0; i<psize(); i++)
	{
		if (ip[i]) npar++;
	}
	
	fout << " Number of constraints        : " << setw(5) << varsize() << endl
		 << " Number of refined parameters : " << setw(5) << npar << endl
		 << " Number of fixed parameters   : " << setw(5) << psize() - npar << endl;

	if (npar != 0)
	{
		int i, j;
		
		fout << endl << " Refinement parameters :\n";

		for (i=0, j=0; i<psize(); i++)
		{
			if (ip[i])
				fout << setw(4) << id[i] << ": " << cc(p[i],dp[i]);
			else
				fout << setw(4) << id[i] << ": " << cc(p[i],0);

			j++;
			if (j%3) fout << "  ";
			else fout << endl;
		}
		if (j%3) fout << endl;

		fout << " " << string(78,'-') << endl 
			 << " REFINEMENT INFORMATION :" << endl
			 << " " << string(78,'-') << endl;

		fout << " Number of iterations : " << iter << endl;

		fout << " Reduced chi squared    : " << setw(15) << redchisq << endl
		     << " Rw - value             : " << setw(15) << fit_rw << endl;

		fout << endl << " Correlations greater than 0.8 : \n\n";

		bool lkor = false;

		for (i=0; i<psize(); i++)
		{
			if (!ip[i]) continue;

			for (j=i+1; j<psize(); j++)
			{
				if (!ip[j]) continue;
				
				double corr = covar[i][j]/dp[i]/dp[j];
				
				if (abs(corr) > 0.8)
				{
					fout << "   Corr(p[" << id[i] << "], p[" << id[j] << "]) = " << corr << endl;
					lkor = true;
				}
			}
		}
		if (!lkor) 
			fout << "   *** none ***\n";
	}
}
