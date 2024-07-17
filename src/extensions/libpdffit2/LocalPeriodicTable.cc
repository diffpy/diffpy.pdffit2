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
*           factors.  This class throws ValueError for unknown elements.
*
***********************************************************************/

#include <stdexcept>
#include "LocalPeriodicTable.h"
#include "exceptions.h"

using namespace std;


// public class methods


const LocalPeriodicTable* LocalPeriodicTable::instance()
{
    static unique_ptr<LocalPeriodicTable> the_table(new LocalPeriodicTable());
    return the_table.get();
}


// constructor


LocalPeriodicTable::LocalPeriodicTable()
{
    // private data
    this->_periodic_table = PeriodicTable::instance();
}


// public methods


const AtomType* LocalPeriodicTable::name(const string& nm) const
{
    const AtomType* atp0;
    try {
        atp0 = this->_periodic_table->name(nm);
    }
    catch (runtime_error e) {
        throw ValueError(e.what());
    }
    const AtomType* atp = local_symbol(atp0->symbol);
    return atp;
}


const AtomType* LocalPeriodicTable::symbol(const string& smbl) const
{
    const AtomType* atp = local_symbol(smbl);
    return atp;
}


const AtomType* LocalPeriodicTable::lookup(const string& pat) const
{
    const AtomType* atp0;
    try {
        atp0 = this->_periodic_table->lookup(pat);
    }
    catch (runtime_error e) {
        throw ValueError(e.what());
    }
    const AtomType* atp = local_symbol(atp0->symbol);
    return atp;
}


bool LocalPeriodicTable::has(const string& pat) const
{
    return this->_periodic_table->has(pat);
}


void LocalPeriodicTable::reset(const string& smbl)
{
    const AtomType* atp0;
    try {
        atp0 = this->_periodic_table->symbol(smbl);
    }
    catch (runtime_error e) {
        throw ValueError(e.what());
    }
    // overwrite with default only when smbl is in the local table
    if (this->_local_table.count(smbl))
    {
        this->_local_table[smbl] = *atp0;
    }
}


void LocalPeriodicTable::setXsf(const string& smbl, double xsf)
{
    AtomType* atp = this->local_symbol(smbl);
    atp->xsf = xsf;
}


void LocalPeriodicTable::setNsf(const string& smbl, double nsf)
{
    AtomType* atp = this->local_symbol(smbl);
    atp->nsf = nsf;
}


// private methods


AtomType* LocalPeriodicTable::local_symbol(const string& smbl) const
{
    map<string,AtomType>::iterator atplocal;
    atplocal = this->_local_table.find(smbl);
    if (atplocal == this->_local_table.end())
    {
        const AtomType* atp0;
        try {
            atp0 = this->_periodic_table->symbol(smbl);
        }
        catch (runtime_error e) {
            throw ValueError(e.what());
        }
        atplocal = this->_local_table.insert(make_pair(smbl, *atp0)).first;
    }
    AtomType* atp = &(atplocal->second);
    return atp;
}


// End of file
