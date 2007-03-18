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
* Output methods for Phase, DataSet and Fit classes
*
* Comments:
*
* $Id$
*
***********************************************************************/

// ensure math constants get defined for MSVC
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "pdffit.h"
#include "StringUtils.h"

/************************************************************
 *  Outputs phase information for phase 'ipha' on file 'id'
 *  Thu Oct 13 2005 - CLF Modified code to handle
 *  general types of streams.
 *************************************************************/
void Phase::output(ostream &fout)
{
    //test();
    fout << " " << string(78,'-') << endl
	<< " PHASE " << iphase << " : " << name << endl
	<< " " << string(78,'-') << endl;

    fout << " Scale factor          : " << cc(skal,dskal) << endl;

    if (corr_max > 0.0)
	fout << " Correlation limit [A] : " << corr_max << endl;

    fout << " Quad. corr. factor    : " << cc(delta2,ddelta2) << endl;

    fout << " Lin. corr. factor     : " << cc(delta1,ddelta1) << endl;

    fout << " Low r sigma ratio     : " << cc(srat,dsrat) << endl
	<< " R cutoff [A]          : " << cc(rcut,0.0) << endl;

    fout << " Lattice parameters    : ";
    int i;
    for (i=0; i<3; i++)
	fout << cc(a0[i], da0[i]);
    fout << endl;

    fout << "           & angles    : ";
    for (i=0; i<3; i++)
	fout << cc(win[i], dwin[i]);
    fout << endl << endl;

    fout << " Atom positions & occupancies : " << endl;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	fout << "   " << setw(4) << toupper(ai->atom_type->symbol);
	for (i=0; i<3; i++)	fout << cc(ai->pos[i],ai->dpos[i]);
	fout << cc(ai->occ,ai->docc) << endl;
    }
    fout << endl;

    fout << " Anisotropic temperature factors : " << endl;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	fout << "   " << setw(4) << toupper(ai->atom_type->symbol);
	for (i=0; i<3; i++)	fout << cc(ai->u[i],ai->du[i]);
	fout << endl;
	if ( ai->u[3] || ai->u[4] || ai->u[5])
	{
	    for (i=3; i<6; i++)
	    {
		fout << "            " << cc(ai->u[i],ai->du[i]);
	    }
	    fout << endl;
	}
    }
}

string DataSet::selectedAtomsString(int ip, char ijchar)
{
    if (!psel[ip])  return string("");
    if (ijchar != 'i' && ijchar != 'j')
    {
	ostringstream emsg;
	emsg << "Invalid value of ijchar '" << ijchar << "'";
	throw ValueError(emsg.str());
    }
    // build string of selected indices per each atom type
    // also check if any type is selected and ignored at the same time
    map<AtomType*, string> selidxstr;
    set<AtomType*> ignored_types;
    Phase* ph = psel[ip];
    set<int>& ignored = ijchar == 'i' ? phase_ignore_i[ph] : phase_ignore_j[ph];
    for (int aidx = 0; aidx < ph->natoms; ++aidx)
    {
	AtomType* atp = ph->atom[aidx].atom_type;
	if (ignored.count(aidx))	ignored_types.insert(atp);
	else
	{
	    ostringstream sidx;
	    sidx << ' ' << aidx + 1;
	    selidxstr[atp] += sidx.str();
	}
    }
    ostringstream ssel;
    for (   vector<AtomType*>::iterator atp = ph->atom_types.begin();
	    atp != ph->atom_types.end(); ++atp )
    {
	if (!selidxstr.count(*atp))	continue;
	ssel << "  " << toupper((*atp)->symbol);
	if (ignored_types.count(*atp))	ssel << selidxstr[*atp];
    }
    return ssel.str();
}

/***********************************************************
  Outputs information about this data set 
  Thu Oct 13 2005 - CLF
  Modified code to handle general types of
  streams.
************************************************************/
void DataSet::output(ostream& fout)
{
    fout << " " << string(78,'-') << endl
	<< " DATA SET : " << iset << " (" << name << ")" << endl
	<< " " << string(78,'-') << endl;

    fout << " Data range in r [A]   : " << setw(8) << rmin << " -> " << setw(8) << rmax
	<< "      Step dr  : " << setw(8) << deltar << endl;

    fout << " Calculated range      : " << setw(8) << rcmin << " -> " << setw(8) << rcmax << endl;

    fout << " Refinement r range    : " << setw(8) << rfmin << " -> " << setw(8) << rfmax 
	<< "      Data pts : " << setw(5) << nfmin << " -> " << setw(5) << nfmax << endl;

    fout << endl << " Experimental settings : " << endl;

    if (scattering_type == 'X')
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

    for(size_t ip = 0; ip != psel.size(); ip++)
    {
	if (!psel[ip])	continue;
	fout << "   Phase " << ip+1 << " :\n";

	fout << "     Atoms (i) :";
	fout << selectedAtomsString(ip, 'i') << endl;
	fout << "     Atoms (j) :";
	fout << selectedAtomsString(ip, 'j') << endl;
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

		if (fabs(corr) > 0.8)
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

// End of file
