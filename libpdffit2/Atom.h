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
* $Id$
*
***********************************************************************/

#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED

#include <iostream>

class AtomType;

class PdfFit;
class DataSet;
class Phase;

class Atom
{
    private:

	enum AtomFormat { DISCUS, PDFFIT };
	static AtomFormat streamformat;
	std::istream& read_discus_atom(std::istream& in);
	std::istream& read_pdffit_atom(std::istream& in);
	int offset;
	// and guys who use offset
	friend class PdfFit;
	friend class DataSet;
	friend class Phase;

    // as atoms data are often accessed by other functions we make its data public
    public:

	AtomType* atom_type;
	double weight;	// normalized scattering factor
	double pos[3], dpos[3];
	double u[6], du[6];
	double occ, docc;

	void setDiscusFormat() { streamformat = DISCUS; }
	void setPdffitFormat() { streamformat = PDFFIT; }
	friend std::istream& operator>>(std::istream& in, Atom& a);

};

std::istream& operator>>(std::istream& in, Atom& a);

#endif	// ATOM_H_INCLUDED
