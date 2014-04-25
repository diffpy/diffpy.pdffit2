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
* class AtomType
*
* Comments: storage of element properties like name, symbol,
*     atomic mass or scattering power
*
***********************************************************************/

#ifndef ATOMTYPE_H_INCLUDED
#define ATOMTYPE_H_INCLUDED

#include <string>
#include <sstream>
#include <stdexcept>

class AtomType
{
    public:
	// Data
	std::string symbol; // element symbol
	std::string name;   // element name
	bool isotope;	    // flag for isotopes
	int z;		    // atomic number
	double M;	    // atomic mass
	int charge;	    // total charge in e
	double radius;	    // ionic radius
	double xsf;	    // x-ray scattering factor
	double nsf;	    // neutron scattering factor

	// Constructor
	AtomType()
	{
	    z = 0;
	    M = 0.0;
	    isotope = false;
	    charge = 0;
	    radius = 0.0;
	    xsf = 0.0;
	    nsf = 0.0;
	}

	// Methods
	double sf(char scattering_type) const
	{
	    switch (scattering_type)
	    {
		case 'x':
		case 'X':
		    return xsf;
		case 'n':
		case 'N':
		    return nsf;
		default:
		    std::ostringstream emsg;
                    emsg << "Invalid scattering type '" <<
                        scattering_type << "'";
		    throw std::runtime_error(emsg.str());
	    }
	    return 0.0;
	}
};

#endif	// ATOMTYPE_H_INCLUDED
