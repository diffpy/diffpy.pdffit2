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
    FormatValueWithStd value_std;
    value_std.left();

    fout << " " << string(78,'-') << '\n'
	<< " PHASE " << iphase << " : " << name << '\n'
	<< " " << string(78,'-') << endl;

    fout << " Scale factor          : " << value_std(skal, dskal) << endl;

    if (corr_max > 0.0)
    {
	fout << " Correlation limit [A] : " << corr_max << endl;
    }

    fout << " Quad. corr. factor    : "
	<< value_std(delta2, ddelta2) << endl;

    fout << " Lin. corr. factor     : "
	<< value_std(delta1, ddelta1) << endl;

    fout << " Low r sigma ratio     : " << value_std(sratio, dsratio) << '\n'
	<<  " R cutoff [A]          : " << value_std(rcut, 0.0) << endl;

    value_std.leading_blank(true).left();
    fout << " Lattice parameters    :"
	<< value_std.width(20)(a0[0], da0[0])
	<< value_std.width(20)(a0[1], da0[1])
	<< value_std.width(0)(a0[2], da0[2]) << '\n';

    fout << "           & angles    :"
	<< value_std.width(20)(win[0], dwin[0])
	<< value_std.width(20)(win[1], dwin[1])
	<< value_std.width(0)(win[2], dwin[2]) << '\n';

    fout << endl;

    fout << " Atom positions & occupancies :" << endl;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	fout << "   "
	    << setw(4) << left << toupper(ai->atom_type->symbol) << setw(0);
	value_std.width(20);
	for (int i = 0; i < 3; i++)
	{
	    fout << value_std(ai->pos[i], ai->dpos[i]);
	}
	fout << value_std.width(0)(ai->occ,ai->docc) << endl;
    }
    fout << endl;

    fout << " Anisotropic temperature factors :" << endl;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	fout << "   " << setw(4) << toupper(ai->atom_type->symbol);
	for (int i = 0; i < 3; i++)
	{
	    fout << value_std.width(i<2 ? 20 : 0)(ai->u[i],ai->du[i]);
	}
	fout << endl;
	if ( ai->u[3] || ai->u[4] || ai->u[5])
	{
	    for (int i = 0; i < 3; i++)
	    {
		fout << "            ";
		fout << value_std.width(i<2 ? 20 : 0)(ai->u[i],ai->du[i]);
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
    FormatValueWithStd value_std;

    fout << " " << string(78,'-') << '\n'
	<< " DATA SET : " << iset << " (" << name << ")" << '\n'
	<< " " << string(78,'-') << endl;

    fout << " Data range in r [A]   : " << setw(8) << left << rmin << " -> "
	<< setw(8) << rmax << "      Step dr  : " << setw(0) << deltar << endl;

    fout << " Calculated range      : " << setw(8) << rcmin << " -> "
	<< setw(0) << rcmax << endl;

    fout << " Refinement r range    : " << setw(8) << rfmin << " -> "
	<< setw(8) << rfmax
	<< "      Data pts : " << setw(5) << nfmin << " -> "
	<< setw(0) << nfmax << endl;

    fout << endl << " Experimental settings :" << endl;

    if (scattering_type == 'X')
	fout << "   Radiation           : X-Rays\n";
    else
	fout << "   Radiation           : Neutrons\n";

    if (qmax <= 0.0)
	fout << "   Termination at Qmax : not applied\n";
    else
	fout << "   Termination at Qmax : " << qmax << " A**-1\n";

    if (qdamp <= 0.0)
	fout << "   DQ dampening Qdamp  : not applied\n";
    else
	fout << "   DQ dampening Qdamp  : " << value_std(qdamp, dqdamp) << " A**-1\n";

    if (qbroad <= 0.0)
	fout << "   DQ broadening Qbroad: not applied\n";
    else
	fout << "   DQ broadening Qbroad: " << value_std(qbroad, dqbroad) << " A**-1\n";

    if (spdiameter <= 0.0)
	fout << "   Particle diameter   : not applied\n";
    else
	fout << "   Particle diameter   : " << value_std(spdiameter, dspdiameter) << " A\n";

    fout << "   Scale factor        : " << value_std(skal, dskal) << endl;

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
    fout << " " << string(78,'-') << '\n'
	<< " PARAMETER INFORMATION :" << '\n'
	<< " " << string(78,'-') << endl;

    int npar = 0;

    for (int i=0; i<psize(); i++)
    {
	if (ip[i]) npar++;
    }

    fout << left << setw(0)
	<< " Number of constraints        : " << varsize() << '\n'
	<< " Number of refined parameters : " << npar << '\n'
	<< " Number of fixed parameters   : " << psize() - npar << endl;

    if (npar != 0)
    {
	fout << endl << " Refinement parameters :\n";

	FormatValueWithStd value_std;
	value_std.leading_blank(true).left();
	for (int i = 0; i < psize(); i++)
	{
	    bool lastword = (i + 1) % 3 == 0 || i + 1 == psize();
	    double dpi = ip[i] ? dp[i] : 0.0;
	    fout << right << ' ' << setw(3) << id[i] << setw(0) << ":"
		<< value_std.width(lastword ? 0 : 20)(p[i], dpi);
	    fout << (lastword ? "\n" : "  ");
	}

	fout << " " << string(78,'-') << '\n'
	    << " REFINEMENT INFORMATION:\n"
	    << " " << string(78,'-') << endl;

	fout << " Number of iterations : " << iter << endl;

	fout << " Reduced chi squared  : " << redchisq << '\n'
	    <<  " Rw - value           : " << fit_rw << '\n' << endl;

	fout << " Correlations greater than 0.8 :\n\n";

	bool lkor = false;

	for (int i = 0; i < psize(); i++)
	{
	    if (!ip[i]) continue;

	    for (int j = i + 1; j < psize(); j++)
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
