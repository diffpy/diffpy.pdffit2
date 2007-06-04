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
* class LocalPeriodicTable
*
* Comments: same as PeriodicTable, but throws ValueError exception
*
* $Id$
*
***********************************************************************/

#ifndef LOCALPERIODICTABLE_H_INCLUDED
#define LOCALPERIODICTABLE_H_INCLUDED

#include <stdexcept>

#include "exceptions.h"
#include "PeriodicTable.h"

class LocalPeriodicTable : public PeriodicTable
{
    public:

	void defAtomType(const AtomType& atp)
	{
	    try {
		PeriodicTable::defAtomType(atp);
	    }
	    catch (std::runtime_error e) {
		throw ValueError(e.what());
	    }
	}

	void reset(AtomType* atp)
	{
	    try {
		PeriodicTable::reset(atp);
	    }
	    catch (std::runtime_error e) {
		throw ValueError(e.what());
	    }
	}

	// not very clean, but seems to be working fine
	static LocalPeriodicTable* instance()
	{
	    PeriodicTable* pt = PeriodicTable::instance();
	    return static_cast<LocalPeriodicTable*>(pt);
	};

};

#endif	// LOCALPERIODICTABLE_H_INCLUDED
