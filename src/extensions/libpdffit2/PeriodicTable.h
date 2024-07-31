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
* class PeriodicTable
*
* Comments: singleton class, use PeriodicTable::instance()
*	    for its pointer
*
***********************************************************************/

#ifndef PERIODICTABLE_H_INCLUDED
#define PERIODICTABLE_H_INCLUDED

#include <deque>
#include <string>
#include <map>
#include <memory>
#include "AtomType.h"
#include "support_msvc9.h"

// PeriodicTable
class PeriodicTable
{
    private:

	// PeriodicTable is a singleton class
	PeriodicTable();

    public:

	// Access to singleton instance
	static PeriodicTable* instance()
	{
	    static std::unique_ptr<PeriodicTable> the_table(new PeriodicTable());
	    return the_table.get();
	}
	// Destructor
	~PeriodicTable();

	// Methods
	AtomType* name(const std::string& s);
	AtomType* symbol(const std::string& s);
	AtomType* lookup(std::string s);    // icase lookup
	bool has(const std::string& s);
	void defAtomType(const AtomType& atp);
	void deleteAtomType(const AtomType* atp);
	void reset(AtomType* atp);	    // retrieve atp from pt_backup
	void resetAll();		    // reset all elements

    private:

	// Data Members
	std::map<std::string,AtomType*> name_index;
	std::map<std::string,AtomType*> symbol_index;
	std::deque<AtomType*> pt_public;
	std::deque<AtomType*> pt_backup;

	// Methods
	void init();
	void clear();
	void fill_pt_backup();
};

#endif	// PERIODICTABLE_H_INCLUDED
