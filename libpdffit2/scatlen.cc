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
* Phase methods for accessing scattering factors.
*
* Comments: Up to date with 1.3.10 Fortran version.  
*	    In Fortran this was fourier.f
*
* $Id$
*
***********************************************************************/

#include <iomanip>
#include <sstream>

#include "LocalPeriodicTable.h"
#include "StringUtils.h"
#include "pdffit.h"

using NS_PDFFIT2::pout;

void Phase::show_scat(char tp)
{
    *pout << get_scat_string(tp);
}

string Phase::get_scat_string(char tp)
{
    stringstream sout;
    vector<AtomType*>::iterator atp;
    for (atp = atom_types.begin(); atp != atom_types.end(); ++atp)
    {
        sout << get_scat_string(tp, *atp);
    }
    return sout.str();
}

string Phase::get_scat_string(char tp, string symbol)
{
    LocalPeriodicTable* pt = LocalPeriodicTable::instance();
    return get_scat_string(tp, pt->lookup(symbol));
}

string Phase::get_scat_string(char tp, AtomType* atp)
{
    stringstream sout;
    string usymbol = toupper(atp->symbol);
    switch(tp)
    {
	case 'N':
	    sout << "Neutron scattering length for " << usymbol << " :  ";
	    sout << atp->nsf << endl;
	    break;
	case 'X':
	    sout << "X-ray scattering factor for " << usymbol << " :  ";
	    sout << atp->xsf << endl;
	    break;
    }
    return sout.str();
}

void Phase::set_scat(char tp, string symbol, double value)
{
    LocalPeriodicTable* pt = LocalPeriodicTable::instance();
    AtomType* atp = pt->lookup(symbol);
    switch (tp)
    {
	case 'X':   atp->xsf = value;
		    break;
	case 'N':   atp->nsf = value;
		    break;
    }
    *pout << get_scat_string(tp, atp);
}

void Phase::reset_scat(char tp, string symbol)
{
    LocalPeriodicTable* pt = LocalPeriodicTable::instance();
    AtomType* atp = pt->lookup(symbol);
    pt->reset(atp);
    *pout << get_scat_string(tp, atp);
}

// End of file
