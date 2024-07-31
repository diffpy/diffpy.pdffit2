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

#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED

#include <iostream>
#include "AtomType.h"

class PdfFit;
class DataSet;
class Phase;

class Atom
{
    // friends who need to touch offset
    friend class PdfFit;
    friend class DataSet;
    friend class Phase;
    friend std::istream& operator>>(std::istream& in, Atom& a);

    public:

        // class methods
	static void setDiscusFormat();
	static void setPdffitFormat();

        // data
	const AtomType* atom_type;
	double weight;	// normalized scattering factor
	double pos[3], dpos[3];
	double u[6], du[6];
	double occ, docc;

    private:

        // types and class data
	enum AtomFormat { DISCUS, PDFFIT };
	static AtomFormat streamformat;

        // data
	int offset;

        // methods
	std::istream& read_discus_atom(std::istream& in);
	std::istream& read_pdffit_atom(std::istream& in);
};

// non-member operators

std::istream& operator>>(std::istream& in, Atom& a);

#endif	// ATOM_H_INCLUDED
