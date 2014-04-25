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
    vector<const AtomType*>::iterator atp;
    for (atp = atom_types.begin(); atp != atom_types.end(); ++atp)
    {
        sout << get_scat_string(tp, *atp);
    }
    return sout.str();
}

string Phase::get_scat_string(char tp, string smbpat)
{
    const LocalPeriodicTable* lpt = getPeriodicTable();
    const AtomType* atp = lpt->lookup(smbpat);
    return get_scat_string(tp, atp);
}

string Phase::get_scat_string(char tp, const AtomType* atp)
{
    stringstream sout;
    string usymbol = toupper(atp->symbol);
    switch(tp)
    {
	case 'n':
	case 'N':
	    sout << "Neutron scattering length for " << usymbol << " :  ";
	    break;
	case 'x':
	case 'X':
	    sout << "X-ray scattering factor for " << usymbol << " :  ";
	    break;
    }
    // this also throws runtime_error for invalid tp value
    sout << atp->sf(tp);
    return sout.str();
}

void Phase::set_scat(char tp, const string& smbpat, double value)
{
    LocalPeriodicTable* lpt = getPeriodicTable();
    const string& stdsmbl = lpt->lookup(smbpat)->symbol;
    switch (tp)
    {
	case 'n':
	case 'N':
            lpt->setNsf(stdsmbl, value);
            break;
	case 'x':
	case 'X':
            lpt->setXsf(stdsmbl, value);
            break;
        default:
            ostringstream emsg;
            emsg << "Invalid scattering type '" << tp << "'";
            throw runtime_error(emsg.str());
    }
    const AtomType* atp = lpt->symbol(stdsmbl);
    *pout << get_scat_string(tp, atp);
}

void Phase::reset_scat(const string& smbpat)
{
    LocalPeriodicTable* lpt = getPeriodicTable();
    const AtomType* atp = lpt->lookup(smbpat);
    const string& stdsmbl = atp->symbol;
    lpt->reset(stdsmbl);
    *pout << get_scat_string('N', stdsmbl);
    *pout << get_scat_string('X', stdsmbl);
}


LocalPeriodicTable* Phase::getPeriodicTable()
{
    return &this->_local_periodic_table;
}


// End of file
