/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2008 trustees of the Michigan State University
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
* Comments: Light-weight class which allows redefinitions of scattering
*           factors.  It keeps a local copy of each looked up species,
*           which can be redefined.
*           This class also throws ValueError for unknown elements.
*
***********************************************************************/

#ifndef LOCALPERIODICTABLE_H_INCLUDED
#define LOCALPERIODICTABLE_H_INCLUDED

#include <map>
#include <string>
#include "PeriodicTable.h"

class LocalPeriodicTable
{
    public:

        // class methods

        // common instance for clients that do not need to own one
	static const LocalPeriodicTable* instance();

        // constructor
        LocalPeriodicTable();

        // methods
	const AtomType* name(const std::string& nm) const;
	const AtomType* symbol(const std::string& smbl) const;
	const AtomType* lookup(const std::string& pat) const;
	bool has(const std::string& pat) const;
	void reset(const std::string& smbl);
        void setXsf(const std::string& smbl, double xsf);
        void setNsf(const std::string& smbl, double nsf);

    private:

        // data
        PeriodicTable* _periodic_table;
	mutable std::map<std::string,AtomType> _local_table;

        // methods
        AtomType* local_symbol(const std::string& smbl) const;

};

#endif	// LOCALPERIODICTABLE_H_INCLUDED
