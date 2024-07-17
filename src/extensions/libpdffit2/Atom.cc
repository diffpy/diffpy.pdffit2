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
* class Atom
*
* Comments:
*
***********************************************************************/

#include <cmath>

#include "Atom.h"
#include "MathUtils.h"
#include "LocalPeriodicTable.h"

using namespace std;

// class data - private


Atom::AtomFormat Atom::streamformat = Atom::DISCUS;


// public class methods


void Atom::setDiscusFormat()
{
    Atom::streamformat = DISCUS;
}


void Atom::setPdffitFormat()
{
    Atom::streamformat = PDFFIT;
}


// private methods


istream& Atom::read_discus_atom(istream& in)
{
    double B;
    string symbol;
    const double fac = 1.0 / (8.0*M_PI*M_PI);
    in >> symbol >> pos[0] >> pos[1] >> pos[2] >> B;
    if (!in)	return in;
    // here we read successfully
    const LocalPeriodicTable* lpt = LocalPeriodicTable::instance();
    atom_type = lpt->lookup(symbol);
    fill_n(u, 3, fac*B);
    fill_n(u+3, 3, 0.0);
    occ = 1.0;
    docc = 0.0;
    fill_n(dpos, 3, 0.0);
    fill_n(du, 6, 0.0);
    return in;
}


istream& Atom::read_pdffit_atom(istream& in)
{
    string symbol;
    in >> symbol >> pos[0] >> pos[1] >> pos[2] >> occ >>
	dpos[0] >> dpos[1] >> dpos[2] >> docc >>
	u[0] >> u[1] >> u[2] >>
	du[0] >> du[1] >> du[2] >>
	u[3] >> u[4] >> u[5] >>
	du[3] >> du[4] >> du[5];
    if (!in)	return in;
    // here we read successfully
    const LocalPeriodicTable* lpt = LocalPeriodicTable::instance();
    atom_type = lpt->lookup(symbol);
    return in;
}


// non-member operators


istream& operator>>(istream& in, Atom& a)
{
    switch (Atom::streamformat)
    {
	case Atom::DISCUS:	return a.read_discus_atom(in);
	case Atom::PDFFIT:	return a.read_pdffit_atom(in);
    };
    return in;
}


// End of file
