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

#include <sstream>
#include <stdexcept>
#include <locale>

// MS compatibility fix
#include <algorithm>

#include "PeriodicTable.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// helper function
////////////////////////////////////////////////////////////////////////
namespace {

inline bool isalphachar(const char& c)
{
    return isalpha(c);
}

} // namespace

////////////////////////////////////////////////////////////////////////
// PeriodicTable definitions
////////////////////////////////////////////////////////////////////////

PeriodicTable::PeriodicTable()
{
    init();
}

PeriodicTable::~PeriodicTable()
{
    clear();
}

AtomType* PeriodicTable::name(const string& s)
{
    map<string,AtomType*>::iterator ii;
    ii = name_index.find(s);
    if (ii == name_index.end())
    {
	ostringstream emsg;
	emsg << "Element or isotope '" << s << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    return ii->second;
}

AtomType* PeriodicTable::symbol(const string& s)
{
    map<string,AtomType*>::iterator ii;
    ii = symbol_index.find(s);
    if (ii == symbol_index.end())
    {
	ostringstream emsg;
	emsg << "Element or isotope '" << s << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    return ii->second;
}

AtomType* PeriodicTable::lookup(string s)
{
    // force standard case
    string::iterator sii = find_if(s.begin(), s.end(), isalphachar);
    if (sii != s.end())
    {
	*sii = toupper(*sii);
	for (sii++; sii != s.end(); ++sii)  *sii = tolower(*sii);
    }
    map<string,AtomType*>::iterator ii;
    ii = symbol_index.find(s);
    if (    ii == symbol_index.end() &&
	    (ii = name_index.find(s)) == name_index.end() )
    {
	ostringstream emsg;
	emsg << "Element  or isotope '" << s << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    return ii->second;
}

bool PeriodicTable::has(const std::string& s)
{
    bool rv;
    try
    {
	lookup(s);
	rv = true;
    }
    catch(runtime_error)
    {
	rv = false;
    }
    return rv;
}

void PeriodicTable::defAtomType(const AtomType& atp)
{
    // check if already defined
    if (symbol_index.count(atp.symbol))
    {
	ostringstream emsg;
	emsg << "Element symbol '" << atp.symbol
	    << "' already defined.";
	throw runtime_error(emsg.str());
    }
    if (name_index.count(atp.name))
    {
	ostringstream emsg;
	emsg << "Element name '" << atp.name
	    << "' already defined.";
	throw runtime_error(emsg.str());
    }
    // all should be fine here
    pt_backup.push_back(new AtomType(atp));
    AtomType* newatp = new AtomType(atp);
    pt_public.push_back(newatp);
    symbol_index[newatp->symbol] = newatp;
    name_index[newatp->name] = newatp;
}

void PeriodicTable::deleteAtomType(const AtomType* atp)
{
    deque<AtomType*>::iterator ii;
    ii = find(pt_public.begin(), pt_public.end(), atp);
    if (ii == pt_public.end())	return;
    // here we need to free and remove atp related data
    symbol_index.erase(atp->symbol);
    name_index.erase(atp->name);
    size_t idx = ii - pt_public.begin();
    delete pt_public[idx];  pt_public.erase(pt_public.begin() + idx);
    delete pt_backup[idx];  pt_backup.erase(pt_backup.begin() + idx);
}

void PeriodicTable::reset(AtomType* atp)
{
    deque<AtomType*>::iterator ii;
    if (!count(pt_backup.begin(), pt_backup.end(), atp))
    {
	ostringstream emsg;
	emsg << "Element '" << atp->symbol << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    size_t idx = ii - pt_backup.begin();
    *(pt_public[idx]) = *(pt_backup[idx]);
}

void PeriodicTable::resetAll()
{
    deque<AtomType*>::iterator iipub = pt_public.begin();
    deque<AtomType*>::iterator iibak = pt_backup.begin();
    for ( ; iipub != pt_public.end(); ++iipub, ++iibak)
    {
	**iipub = **iibak;
    }
}

void PeriodicTable::init()
{
    // load element data to pt_backup
    fill_pt_backup();
    // make public copy and initialize lookup maps
    pt_public.resize(pt_backup.size());
    // fill in lookup map
    deque<AtomType*>::iterator iipub = pt_public.begin();
    deque<AtomType*>::iterator iibak = pt_backup.begin();
    for ( ; iipub != pt_public.end(); ++iipub, ++iibak)
    {
	*iipub = new AtomType(**iibak);
	symbol_index[(*iipub)->symbol] = *iipub;
	name_index[(*iipub)->name] = *iipub;
    }
    // add standard symbols for deuterium and tritium
    symbol_index["2-H"] = lookup("D");
    symbol_index["3-H"] = lookup("T");
}

void PeriodicTable::clear()
{
    deque<AtomType*>::iterator iipub = pt_public.begin();
    deque<AtomType*>::iterator iibak = pt_backup.begin();
    for ( ; iipub != pt_public.end(); ++iipub, ++iibak)
    {
	delete *iipub; *iipub = NULL;
	delete *iibak; *iibak = NULL;
    }
    pt_public.clear();
    pt_backup.clear();
}

void PeriodicTable::fill_pt_backup()
{
    // Refs:
    // 1. Albert-Jose Dianoux, Gerry Lander, Neutron Data Booklet,
    //    Second Edition, ILL 2003
    // 2. ionic radii: http://www.fhi-berlin.mpg.de/th/balsac/balm.47.html
    AtomType* atp;
    AtomType* itp;
    // hydrogen
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "H";
	atp->name = "Hydrogen";
	atp->z = 1;
	atp->M = 1.007947;
	atp->radius = 0.4350;
	atp->xsf = 1.0;
	atp->nsf = -3.7409;
	// isotopes
	// 1-hydrogen
	itp = new AtomType(*atp);
	itp->symbol = "1-H";
	itp->name = "1-Hydrogen";
	itp->isotope = true;
	itp->M = 1.0078250321;
	itp->nsf = -3.7423;
	pt_backup.push_back(itp);
	// deuterium
	itp = new AtomType(*atp);
	itp->symbol = "D";
	itp->name = "Deuterium";
	itp->isotope = true;
	itp->M = 2.0141017780;
	itp->nsf = 6.674;
	pt_backup.push_back(itp);
	// tritium
	itp = new AtomType(*atp);
	itp->symbol = "T";
	itp->name = "Tritium";
	itp->isotope = true;
	itp->M = 3.0160492675;
	itp->nsf = 4.792;
	pt_backup.push_back(itp);
    }
    // helium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "He";
	atp->name = "Helium";
	atp->z = 2;
	atp->M = 4.0026022;
	atp->radius = 1.4000;
	atp->xsf = 2.0;
	atp->nsf = 3.263;
	// isotopes
	// 3-helium
	itp = new AtomType(*atp);
	itp->symbol = "3-He";
	itp->name = "3-Helium";
	itp->isotope = true;
	itp->M = 3.0160293097;
	itp->nsf = 5.74;
	pt_backup.push_back(itp);
	// 4-helium
	itp = new AtomType(*atp);
	itp->symbol = "4-He";
	itp->name = "4-Helium";
	itp->isotope = true;
	itp->M = 4.0026032497;
	itp->nsf = 3.26;
	pt_backup.push_back(itp);
    }
    // lithium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Li";
	atp->name = "Lithium";
	atp->z = 3;
	atp->M = 6.9412;
	atp->radius = 1.5199;
	atp->xsf = 3.0;
	atp->nsf = -1.903;
	// isotopes
	// 6-lithium
	itp = new AtomType(*atp);
	itp->symbol = "6-Li";
	itp->name = "6-Lithium";
	itp->isotope = true;
	itp->M = 6.0151223;
	itp->nsf = 2.0;
	pt_backup.push_back(itp);
	// 7-lithium
	itp = new AtomType(*atp);
	itp->symbol = "7-Li";
	itp->name = "7-Lithium";
	itp->isotope = true;
	itp->M = 7.0160040;
	itp->nsf = -2.22;
	pt_backup.push_back(itp);
    }
    // beryllium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Be";
	atp->name = "Beryllium";
	atp->z = 4;
	atp->M = 9.0121823;
	atp->radius = 1.1430;
	atp->xsf = 4.0;
	atp->nsf = 7.791;
    }
    // boron
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "B";
	atp->name = "Boron";
	atp->z = 5;
	atp->M = 10.8117;
	atp->radius = 0.9750;
	atp->xsf = 5.0;
	atp->nsf = 5.304;
	// isotopes
	// 10-boron
	itp = new AtomType(*atp);
	itp->symbol = "10-B";
	itp->name = "10-Boron";
	itp->isotope = true;
	itp->M = 10.0129370;
	itp->nsf = -0.2;
	pt_backup.push_back(itp);
	// 11-boron
	itp = new AtomType(*atp);
	itp->symbol = "11-B";
	itp->name = "11-Boron";
	itp->isotope = true;
	itp->M = 11.0093055;
	itp->nsf = 6.65;
	pt_backup.push_back(itp);
    }
    // carbon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "C";
	atp->name = "Carbon";
	atp->z = 6;
	atp->M = 12.01078;
	atp->radius = 0.6550;
	atp->xsf = 6.0;
	atp->nsf = 6.648413;
	// isotopes
	// 12-carbon
	itp = new AtomType(*atp);
	itp->symbol = "12-C";
	itp->name = "12-Carbon";
	itp->isotope = true;
	itp->M = 12.0;
	itp->nsf = 6.6535;
	pt_backup.push_back(itp);
	// 13-carbon
	itp = new AtomType(*atp);
	itp->symbol = "13-C";
	itp->name = "13-Carbon";
	itp->isotope = true;
	itp->M = 13.0033548378;
	itp->nsf = 6.19;
	pt_backup.push_back(itp);
    }
    // nitrogen
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "N";
	atp->name = "Nitrogen";
	atp->z = 7;
	atp->M = 14.00672;
	atp->radius = 0.7500;
	atp->xsf = 7.0;
	atp->nsf = 9.36;
	// isotopes
	// 14-nitrogen
	itp = new AtomType(*atp);
	itp->symbol = "14-N";
	itp->name = "14-Nitrogen";
	itp->isotope = true;
	itp->M = 14.0030740052;
	itp->nsf = 9.37;
	pt_backup.push_back(itp);
	// 15-nitrogen
	itp = new AtomType(*atp);
	itp->symbol = "15-N";
	itp->name = "15-Nitrogen";
	itp->isotope = true;
	itp->M = 15.0001088984;
	itp->nsf = 6.44;
	pt_backup.push_back(itp);
    }
    // oxygen
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "O";
	atp->name = "Oxygen";
	atp->z = 8;
	atp->M = 15.99943;
	atp->radius = 0.7300;
	atp->xsf = 8.0;
	atp->nsf = 5.8054;
	// isotopes
	// 16-oxygen
	itp = new AtomType(*atp);
	itp->symbol = "16-O";
	itp->name = "16-Oxygen";
	itp->isotope = true;
	itp->M = 15.9949146221;
	itp->nsf = 5.805;
	pt_backup.push_back(itp);
	// 17-oxygen
	itp = new AtomType(*atp);
	itp->symbol = "17-O";
	itp->name = "17-Oxygen";
	itp->isotope = true;
	itp->M = 16.99913150;
	itp->nsf = 5.6;
	pt_backup.push_back(itp);
	// 18-oxygen
	itp = new AtomType(*atp);
	itp->symbol = "18-O";
	itp->name = "18-Oxygen";
	itp->isotope = true;
	itp->M = 17.9991604;
	itp->nsf = 5.84;
	pt_backup.push_back(itp);
    }
    // fluorine
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "F";
	atp->name = "Fluorine";
	atp->z = 9;
	atp->M = 18.99840325;
	atp->radius = 0.7200;
	atp->xsf = 9.0;
	atp->nsf = 5.65412;
    }
    // neon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ne";
	atp->name = "Neon";
	atp->z = 10;
	atp->M = 20.17976;
	atp->radius = 1.6000;
	atp->xsf = 10.0;
	atp->nsf = 4.5666;
	// isotopes
	// 20-neon
	itp = new AtomType(*atp);
	itp->symbol = "20-Ne";
	itp->name = "20-Neon";
	itp->isotope = true;
	itp->M = 19.9924401759;
	itp->nsf = 4.631;
	pt_backup.push_back(itp);
	// 21-neon
	itp = new AtomType(*atp);
	itp->symbol = "21-Ne";
	itp->name = "21-Neon";
	itp->isotope = true;
	itp->M = 20.99384674;
	itp->nsf = 6.66;
	pt_backup.push_back(itp);
	// 22-neon
	itp = new AtomType(*atp);
	itp->symbol = "22-Ne";
	itp->name = "22-Neon";
	itp->isotope = true;
	itp->M = 21.99138551;
	itp->nsf = 3.87;
	pt_backup.push_back(itp);
    }
    // sodium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Na";
	atp->name = "Sodium";
	atp->z = 11;
	atp->M = 22.9897702;
	atp->radius = 1.8579;
	atp->xsf = 11.0;
	atp->nsf = 3.632;
    }
    // magnesium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Mg";
	atp->name = "Magnesium";
	atp->z = 12;
	atp->M = 24.30506;
	atp->radius = 1.6047;
	atp->xsf = 12.0;
	atp->nsf = 5.3754;
	// isotopes
	// 24-magnesium
	itp = new AtomType(*atp);
	itp->symbol = "24-Mg";
	itp->name = "24-Magnesium";
	itp->isotope = true;
	itp->M = 23.98504190;
	itp->nsf = 5.49;
	pt_backup.push_back(itp);
	// 25-magnesium
	itp = new AtomType(*atp);
	itp->symbol = "25-Mg";
	itp->name = "25-Magnesium";
	itp->isotope = true;
	itp->M = 24.98583702;
	itp->nsf = 3.62;
	pt_backup.push_back(itp);
	// 26-magnesium
	itp = new AtomType(*atp);
	itp->symbol = "26-Mg";
	itp->name = "26-Magnesium";
	itp->isotope = true;
	itp->M = 25.98259304;
	itp->nsf = 4.89;
	pt_backup.push_back(itp);
    }
    // aluminium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Al";
	atp->name = "Aluminium";
	atp->z = 13;
	atp->M = 26.9815382;
	atp->radius = 1.4318;
	atp->xsf = 13.0;
	atp->nsf = 3.4495;
    }
    // silicon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Si";
	atp->name = "Silicon";
	atp->z = 14;
	atp->M = 28.08553;
	atp->radius = 1.1758;
	atp->xsf = 14.0;
	atp->nsf = 4.15071;
	// isotopes
	// 28-silicon
	itp = new AtomType(*atp);
	itp->symbol = "28-Si";
	itp->name = "28-Silicon";
	itp->isotope = true;
	itp->M = 27.9769265327;
	itp->nsf = 4.106;
	pt_backup.push_back(itp);
	// 29-silicon
	itp = new AtomType(*atp);
	itp->symbol = "29-Si";
	itp->name = "29-Silicon";
	itp->isotope = true;
	itp->M = 28.97649472;
	itp->nsf = 4.7;
	pt_backup.push_back(itp);
	// 30-silicon
	itp = new AtomType(*atp);
	itp->symbol = "30-Si";
	itp->name = "30-Silicon";
	itp->isotope = true;
	itp->M = 29.97377022;
	itp->nsf = 4.58;
	pt_backup.push_back(itp);
    }
    // phosphorus
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "P";
	atp->name = "Phosphorus";
	atp->z = 15;
	atp->M = 30.9737612;
	atp->radius = 1.0600;
	atp->xsf = 15.0;
	atp->nsf = 5.131;
    }
    // sulfur
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "S";
	atp->name = "Sulfur";
	atp->z = 16;
	atp->M = 32.0655;
	atp->radius = 1.0200;
	atp->xsf = 16.0;
	atp->nsf = 2.8471;
	// isotopes
	// 32-sulfur
	itp = new AtomType(*atp);
	itp->symbol = "32-S";
	itp->name = "32-Sulfur";
	itp->isotope = true;
	itp->M = 31.97207069;
	itp->nsf = 2.804;
	pt_backup.push_back(itp);
	// 33-sulfur
	itp = new AtomType(*atp);
	itp->symbol = "33-S";
	itp->name = "33-Sulfur";
	itp->isotope = true;
	itp->M = 32.97145850;
	itp->nsf = 4.74;
	pt_backup.push_back(itp);
	// 34-sulfur
	itp = new AtomType(*atp);
	itp->symbol = "34-S";
	itp->name = "34-Sulfur";
	itp->isotope = true;
	itp->M = 33.96786683;
	itp->nsf = 3.48;
	pt_backup.push_back(itp);
	// 36-sulfur
	itp = new AtomType(*atp);
	itp->symbol = "36-S";
	itp->name = "36-Sulfur";
	itp->isotope = true;
	itp->M = 35.96708088;
	itp->nsf = 3.0;
	pt_backup.push_back(itp);
    }
    // chlorine
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cl";
	atp->name = "Chlorine";
	atp->z = 17;
	atp->M = 35.4532;
	atp->radius = 0.9900;
	atp->xsf = 17.0;
	atp->nsf = 9.57928;
	// isotopes
	// 35-chlorine
	itp = new AtomType(*atp);
	itp->symbol = "35-Cl";
	itp->name = "35-Chlorine";
	itp->isotope = true;
	itp->M = 34.96885271;
	itp->nsf = 11.7;
	pt_backup.push_back(itp);
	// 37-chlorine
	itp = new AtomType(*atp);
	itp->symbol = "37-Cl";
	itp->name = "37-Chlorine";
	itp->isotope = true;
	itp->M = 36.96590260;
	itp->nsf = 3.08;
	pt_backup.push_back(itp);
    }
    // argon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ar";
	atp->name = "Argon";
	atp->z = 18;
	atp->M = 39.9481;
	atp->radius = 1.9000;
	atp->xsf = 18.0;
	atp->nsf = 1.9096;
	// isotopes
	// 36-argon
	itp = new AtomType(*atp);
	itp->symbol = "36-Ar";
	itp->name = "36-Argon";
	itp->isotope = true;
	itp->M = 35.96754628;
	itp->nsf = 24.9;
	pt_backup.push_back(itp);
	// 38-argon
	itp = new AtomType(*atp);
	itp->symbol = "38-Ar";
	itp->name = "38-Argon";
	itp->isotope = true;
	itp->M = 37.9627322;
	itp->nsf = 3.5;
	pt_backup.push_back(itp);
	// 40-argon
	itp = new AtomType(*atp);
	itp->symbol = "40-Ar";
	itp->name = "40-Argon";
	itp->isotope = true;
	itp->M = 39.962383123;
	itp->nsf = 1.7;
	pt_backup.push_back(itp);
    }
    // potassium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "K";
	atp->name = "Potassium";
	atp->z = 19;
	atp->M = 39.09831;
	atp->radius = 2.2620;
	atp->xsf = 19.0;
	atp->nsf = 3.672;
	// isotopes
	// 39-potassium
	itp = new AtomType(*atp);
	itp->symbol = "39-K";
	itp->name = "39-Potassium";
	itp->isotope = true;
	itp->M = 38.9637069;
	itp->nsf = 3.79;
	pt_backup.push_back(itp);
	// 40-potassium
	itp = new AtomType(*atp);
	itp->symbol = "40-K";
	itp->name = "40-Potassium";
	itp->isotope = true;
	itp->M = 39.96399867;
	itp->nsf = 3.1;
	pt_backup.push_back(itp);
	// 41-potassium
	itp = new AtomType(*atp);
	itp->symbol = "41-K";
	itp->name = "41-Potassium";
	itp->isotope = true;
	itp->M = 40.96182597;
	itp->nsf = 2.69;
	pt_backup.push_back(itp);
    }
    // calcium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ca";
	atp->name = "Calcium";
	atp->z = 20;
	atp->M = 40.0784;
	atp->radius = 1.9758;
	atp->xsf = 20.0;
	atp->nsf = 4.702;
	// isotopes
	// 40-calcium
	itp = new AtomType(*atp);
	itp->symbol = "40-Ca";
	itp->name = "40-Calcium";
	itp->isotope = true;
	itp->M = 39.9625912;
	itp->nsf = 4.78;
	pt_backup.push_back(itp);
	// 42-calcium
	itp = new AtomType(*atp);
	itp->symbol = "42-Ca";
	itp->name = "42-Calcium";
	itp->isotope = true;
	itp->M = 41.9586183;
	itp->nsf = 3.36;
	pt_backup.push_back(itp);
	// 43-calcium
	itp = new AtomType(*atp);
	itp->symbol = "43-Ca";
	itp->name = "43-Calcium";
	itp->isotope = true;
	itp->M = 42.9587668;
	itp->nsf = -1.56;
	pt_backup.push_back(itp);
	// 44-calcium
	itp = new AtomType(*atp);
	itp->symbol = "44-Ca";
	itp->name = "44-Calcium";
	itp->isotope = true;
	itp->M = 43.9554811;
	itp->nsf = 1.42;
	pt_backup.push_back(itp);
	// 46-calcium
	itp = new AtomType(*atp);
	itp->symbol = "46-Ca";
	itp->name = "46-Calcium";
	itp->isotope = true;
	itp->M = 45.9536928;
	itp->nsf = 3.55;
	pt_backup.push_back(itp);
	// 48-calcium
	itp = new AtomType(*atp);
	itp->symbol = "48-Ca";
	itp->name = "48-Calcium";
	itp->isotope = true;
	itp->M = 47.952534;
	itp->nsf = 0.39;
	pt_backup.push_back(itp);
    }
    // scandium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sc";
	atp->name = "Scandium";
	atp->z = 21;
	atp->M = 44.9559108;
	atp->radius = 1.6545;
	atp->xsf = 21.0;
	atp->nsf = 12.11;
    }
    // titanium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ti";
	atp->name = "Titanium";
	atp->z = 22;
	atp->M = 47.8671;
	atp->radius = 1.4755;
	atp->xsf = 22.0;
	atp->nsf = -3.37013;
	// isotopes
	// 46-titanium
	itp = new AtomType(*atp);
	itp->symbol = "46-Ti";
	itp->name = "46-Titanium";
	itp->isotope = true;
	itp->M = 45.9526295;
	itp->nsf = 4.72;
	pt_backup.push_back(itp);
	// 47-titanium
	itp = new AtomType(*atp);
	itp->symbol = "47-Ti";
	itp->name = "47-Titanium";
	itp->isotope = true;
	itp->M = 46.9517638;
	itp->nsf = 3.53;
	pt_backup.push_back(itp);
	// 48-titanium
	itp = new AtomType(*atp);
	itp->symbol = "48-Ti";
	itp->name = "48-Titanium";
	itp->isotope = true;
	itp->M = 47.9479471;
	itp->nsf = -5.86;
	pt_backup.push_back(itp);
	// 49-titanium
	itp = new AtomType(*atp);
	itp->symbol = "49-Ti";
	itp->name = "49-Titanium";
	itp->isotope = true;
	itp->M = 48.9478708;
	itp->nsf = 0.98;
	pt_backup.push_back(itp);
	// 50-titanium
	itp = new AtomType(*atp);
	itp->symbol = "50-Ti";
	itp->name = "50-Titanium";
	itp->isotope = true;
	itp->M = 49.9447921;
	itp->nsf = 5.88;
	pt_backup.push_back(itp);
    }
    // vanadium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "V";
	atp->name = "Vanadium";
	atp->z = 23;
	atp->M = 50.94151;
	atp->radius = 1.3090;
	atp->xsf = 23.0;
	atp->nsf = -0.44314;
	// isotopes
	// 50-vanadium
	itp = new AtomType(*atp);
	itp->symbol = "50-V";
	itp->name = "50-Vanadium";
	itp->isotope = true;
	itp->M = 49.9471628;
	itp->nsf = 7.6;
	pt_backup.push_back(itp);
	// 51-vanadium
	itp = new AtomType(*atp);
	itp->symbol = "51-V";
	itp->name = "51-Vanadium";
	itp->isotope = true;
	itp->M = 50.9439637;
	itp->nsf = -0.402;
	pt_backup.push_back(itp);
    }
    // chromium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cr";
	atp->name = "Chromium";
	atp->z = 24;
	atp->M = 51.99616;
	atp->radius = 1.2490;
	atp->xsf = 24.0;
	atp->nsf = 3.6357;
	// isotopes
	// 50-chromium
	itp = new AtomType(*atp);
	itp->symbol = "50-Cr";
	itp->name = "50-Chromium";
	itp->isotope = true;
	itp->M = 49.9460496;
	itp->nsf = -4.50;
	pt_backup.push_back(itp);
	// 52-chromium
	itp = new AtomType(*atp);
	itp->symbol = "52-Cr";
	itp->name = "52-Chromium";
	itp->isotope = true;
	itp->M = 51.9405119;
	itp->nsf = 4.914;
	pt_backup.push_back(itp);
	// 0-chromium
	itp = new AtomType(*atp);
	itp->symbol = "53-Cr";
	itp->name = "53-Chromium";
	itp->isotope = true;
	itp->M = 52.9406538;
	itp->nsf = -4.20;
	pt_backup.push_back(itp);
	// 54-chromium
	itp = new AtomType(*atp);
	itp->symbol = "54-Cr";
	itp->name = "54-Chromium";
	itp->isotope = true;
	itp->M = 53.9388849;
	itp->nsf = 4.55;
	pt_backup.push_back(itp);
    }
    // manganese
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Mn";
	atp->name = "Manganese";
	atp->z = 25;
	atp->M = 54.9380499;
	atp->radius = 1.3500;
	atp->xsf = 25.0;
	atp->nsf = -3.75018;
    }
    // iron
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Fe";
	atp->name = "Iron";
	atp->z = 26;
	atp->M = 55.8452;
	atp->radius = 1.2411;
	atp->xsf = 26.0;
	atp->nsf = 9.452;
	// isotopes
	// 54-iron
	itp = new AtomType(*atp);
	itp->symbol = "54-Fe";
	itp->name = "54-Iron";
	itp->isotope = true;
	itp->M = 53.9396148;
	itp->nsf = 4.2;
	pt_backup.push_back(itp);
	// 56-iron
	itp = new AtomType(*atp);
	itp->symbol = "56-Fe";
	itp->name = "56-Iron";
	itp->isotope = true;
	itp->M = 55.9349421;
	itp->nsf = 10.1;
	pt_backup.push_back(itp);
	// 57-iron
	itp = new AtomType(*atp);
	itp->symbol = "57-Fe";
	itp->name = "57-Iron";
	itp->isotope = true;
	itp->M = 56.9353987;
	itp->nsf = 2.3;
	pt_backup.push_back(itp);
	// 58-iron
	itp = new AtomType(*atp);
	itp->symbol = "58-Fe";
	itp->name = "58-Iron";
	itp->isotope = true;
	itp->M = 57.9332805;
	itp->nsf = 15;
	pt_backup.push_back(itp);
    }
    // cobalt
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Co";
	atp->name = "Cobalt";
	atp->z = 27;
	atp->M = 58.9332009;
	atp->radius = 1.2535;
	atp->xsf = 27.0;
	atp->nsf = 2.492;
    }
    // nickel
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ni";
	atp->name = "Nickel";
	atp->z = 28;
	atp->M = 58.69342;
	atp->radius = 1.2460;
	atp->xsf = 28.0;
	atp->nsf = 10.31;
	// isotopes
	// 58-nickel
	itp = new AtomType(*atp);
	itp->symbol = "58-Ni";
	itp->name = "58-Nickel";
	itp->isotope = true;
	itp->M = 57.9353479;
	itp->nsf = 14.4;
	pt_backup.push_back(itp);
	// 60-nickel
	itp = new AtomType(*atp);
	itp->symbol = "60-Ni";
	itp->name = "60-Nickel";
	itp->isotope = true;
	itp->M = 59.9307906;
	itp->nsf = 2.8;
	pt_backup.push_back(itp);
	// 61-nickel
	itp = new AtomType(*atp);
	itp->symbol = "61-Ni";
	itp->name = "61-Nickel";
	itp->isotope = true;
	itp->M = 60.9310604;
	itp->nsf = 7.60;
	pt_backup.push_back(itp);
	// 62-nickel
	itp = new AtomType(*atp);
	itp->symbol = "62-Ni";
	itp->name = "62-Nickel";
	itp->isotope = true;
	itp->M = 61.9283488;
	itp->nsf = -8.7;
	pt_backup.push_back(itp);
	// 64-nickel
	itp = new AtomType(*atp);
	itp->symbol = "64-Ni";
	itp->name = "64-Nickel";
	itp->isotope = true;
	itp->M = 63.9279696;
	itp->nsf = -0.37;
	pt_backup.push_back(itp);
    }
    // copper
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cu";
	atp->name = "Copper";
	atp->z = 29;
	atp->M = 63.5463;
	atp->radius = 1.2780;
	atp->xsf = 29.0;
	atp->nsf = 7.7184;
	// isotopes
	// 63-copper
	itp = new AtomType(*atp);
	itp->symbol = "63-Cu";
	itp->name = "63-Copper";
	itp->isotope = true;
	itp->M = 62.9296011;
	itp->nsf = 6.477;
	pt_backup.push_back(itp);
	// 65-copper
	itp = new AtomType(*atp);
	itp->symbol = "65-Cu";
	itp->name = "65-Copper";
	itp->isotope = true;
	itp->M = 64.9277937;
	itp->nsf = 10.204;
	pt_backup.push_back(itp);
    }
    // zinc
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Zn";
	atp->name = "Zinc";
	atp->z = 30;
	atp->M = 65.4094;
	atp->radius = 1.3325;
	atp->xsf = 30.0;
	atp->nsf = 5.6805;
	// isotopes
	// 64-zinc
	itp = new AtomType(*atp);
	itp->symbol = "64-Zn";
	itp->name = "64-Zinc";
	itp->isotope = true;
	itp->M = 63.9291466;
	itp->nsf = 5.23;
	pt_backup.push_back(itp);
	// 66-zinc
	itp = new AtomType(*atp);
	itp->symbol = "66-Zn";
	itp->name = "66-Zinc";
	itp->isotope = true;
	itp->M = 65.9260368;
	itp->nsf = 5.98;
	pt_backup.push_back(itp);
	// 67-zinc
	itp = new AtomType(*atp);
	itp->symbol = "67-Zn";
	itp->name = "67-Zinc";
	itp->isotope = true;
	itp->M = 66.9271309;
	itp->nsf = 7.58;
	pt_backup.push_back(itp);
	// 68-zinc
	itp = new AtomType(*atp);
	itp->symbol = "68-Zn";
	itp->name = "68-Zinc";
	itp->isotope = true;
	itp->M = 67.9248476;
	itp->nsf = 6.04;
	pt_backup.push_back(itp);
	// 70-zinc
	itp = new AtomType(*atp);
	itp->symbol = "70-Zn";
	itp->name = "70-Zinc";
	itp->isotope = true;
	itp->M = 69.925325;
	itp->nsf = 6.9;
	pt_backup.push_back(itp);
    }
    // gallium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ga";
	atp->name = "Gallium";
	atp->z = 31;
	atp->M = 69.7231;
	atp->radius = 1.3501;
	atp->xsf = 31.0;
	atp->nsf = 7.2882;
	// isotopes
	// 69-gallium
	itp = new AtomType(*atp);
	itp->symbol = "69-Ga";
	itp->name = "69-Gallium";
	itp->isotope = true;
	itp->M = 68.925581;
	itp->nsf = 8.043;
	pt_backup.push_back(itp);
	// 71-gallium
	itp = new AtomType(*atp);
	itp->symbol = "71-Ga";
	itp->name = "71-Gallium";
	itp->isotope = true;
	itp->M = 70.9247050;
	itp->nsf = 6.170;
	pt_backup.push_back(itp);
    }
    // germanium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ge";
	atp->name = "Germanium";
	atp->z = 32;
	atp->M = 72.641;
	atp->radius = 1.2248;
	atp->xsf = 32.0;
	atp->nsf = 8.18520;
	// isotopes
	// 70-germanium
	itp = new AtomType(*atp);
	itp->symbol = "70-Ge";
	itp->name = "70-Germanium";
	itp->isotope = true;
	itp->M = 69.9242504;
	itp->nsf = 10.0;
	pt_backup.push_back(itp);
	// 72-germanium
	itp = new AtomType(*atp);
	itp->symbol = "72-Ge";
	itp->name = "72-Germanium";
	itp->isotope = true;
	itp->M = 71.9220762;
	itp->nsf = 8.51;
	pt_backup.push_back(itp);
	// 73-germanium
	itp = new AtomType(*atp);
	itp->symbol = "73-Ge";
	itp->name = "73-Germanium";
	itp->isotope = true;
	itp->M = 72.9234594;
	itp->nsf = 5.02;
	pt_backup.push_back(itp);
	// 74-germanium
	itp = new AtomType(*atp);
	itp->symbol = "74-Ge";
	itp->name = "74-Germanium";
	itp->isotope = true;
	itp->M = 73.9211782;
	itp->nsf = 7.58;
	pt_backup.push_back(itp);
	// 76-germanium
	itp = new AtomType(*atp);
	itp->symbol = "76-Ge";
	itp->name = "76-Germanium";
	itp->isotope = true;
	itp->M = 75.9214027;
	itp->nsf = 8.2;
	pt_backup.push_back(itp);
    }
    // arsenic
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "As";
	atp->name = "Arsenic";
	atp->z = 33;
	atp->M = 74.921602;
	atp->radius = 1.2000;
	atp->xsf = 33.0;
	atp->nsf = 6.581;
    }
    // selenium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Se";
	atp->name = "Selenium";
	atp->z = 34;
	atp->M = 78.963;
	atp->radius = 1.1600;
	atp->xsf = 34.0;
	atp->nsf = 7.9709;
	// isotopes
	// 74-selenium
	itp = new AtomType(*atp);
	itp->symbol = "74-Se";
	itp->name = "74-Selenium";
	itp->isotope = true;
	itp->M = 73.9224766;
	itp->nsf = 0.8;
	pt_backup.push_back(itp);
	// 76-selenium
	itp = new AtomType(*atp);
	itp->symbol = "76-Se";
	itp->name = "76-Selenium";
	itp->isotope = true;
	itp->M = 75.9192141;
	itp->nsf = 12.2;
	pt_backup.push_back(itp);
	// 77-selenium
	itp = new AtomType(*atp);
	itp->symbol = "77-Se";
	itp->name = "77-Selenium";
	itp->isotope = true;
	itp->M = 76.9199146;
	itp->nsf = 8.25;
	pt_backup.push_back(itp);
	// 78-selenium
	itp = new AtomType(*atp);
	itp->symbol = "78-Se";
	itp->name = "78-Selenium";
	itp->isotope = true;
	itp->M = 77.9173095;
	itp->nsf = 8.24;
	pt_backup.push_back(itp);
	// 80-selenium
	itp = new AtomType(*atp);
	itp->symbol = "80-Se";
	itp->name = "80-Selenium";
	itp->isotope = true;
	itp->M = 79.9165218;
	itp->nsf = 7.48;
	pt_backup.push_back(itp);
	// 82-selenium
	itp = new AtomType(*atp);
	itp->symbol = "82-Se";
	itp->name = "82-Selenium";
	itp->isotope = true;
	itp->M = 81.9167000;
	itp->nsf = 6.34;
	pt_backup.push_back(itp);
    }
    // bromine
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Br";
	atp->name = "Bromine";
	atp->z = 35;
	atp->M = 79.9041;
	atp->radius = 1.1400;
	atp->xsf = 35.0;
	atp->nsf = 6.792;
	// isotopes
	// 79-bromine
	itp = new AtomType(*atp);
	itp->symbol = "79-Br";
	itp->name = "79-Bromine";
	itp->isotope = true;
	itp->M = 78.9183376;
	itp->nsf = 6.79;
	pt_backup.push_back(itp);
	// 81-bromine
	itp = new AtomType(*atp);
	itp->symbol = "81-Br";
	itp->name = "81-Bromine";
	itp->isotope = true;
	itp->M = 80.916291;
	itp->nsf = 6.78;
	pt_backup.push_back(itp);
    }
    // krypton
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Kr";
	atp->name = "Krypton";
	atp->z = 36;
	atp->M = 83.7982;
	atp->radius = 2.0000;
	atp->xsf = 36.0;
	atp->nsf = 7.812;
	// isotopes
	// 86-krypton
	itp = new AtomType(*atp);
	itp->symbol = "86-Kr";
	itp->name = "86-Krypton";
	itp->isotope = true;
	itp->M = 85.9106103;
	itp->nsf = 8.07;
	pt_backup.push_back(itp);
    }
    // rubidium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Rb";
	atp->name = "Rubidium";
	atp->z = 37;
	atp->M = 85.46783;
	atp->radius = 2.4700;
	atp->xsf = 37.0;
	atp->nsf = 7.082;
	// isotopes
	// 85-rubidium
	itp = new AtomType(*atp);
	itp->symbol = "85-Rb";
	itp->name = "85-Rubidium";
	itp->isotope = true;
	itp->M = 84.9117893;
	itp->nsf = 7.07;
	pt_backup.push_back(itp);
	// 87-rubidium
	itp = new AtomType(*atp);
	itp->symbol = "87-Rb";
	itp->name = "87-Rubidium";
	itp->isotope = true;
	itp->M = 86.9091835;
	itp->nsf = 7.27;
	pt_backup.push_back(itp);
    }
    // strontium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sr";
	atp->name = "Strontium";
	atp->z = 38;
	atp->M = 87.621;
	atp->radius = 2.1513;
	atp->xsf = 38.0;
	atp->nsf = 7.022;
	// isotopes
	// 84-strontium
	itp = new AtomType(*atp);
	itp->symbol = "84-Sr";
	itp->name = "84-Strontium";
	itp->isotope = true;
	itp->M = 83.913425;
	itp->nsf = 5.0;
	pt_backup.push_back(itp);
	// 86-strontium
	itp = new AtomType(*atp);
	itp->symbol = "86-Sr";
	itp->name = "86-Strontium";
	itp->isotope = true;
	itp->M = 85.9092624;
	itp->nsf = 5.68;
	pt_backup.push_back(itp);
	// 87-strontium
	itp = new AtomType(*atp);
	itp->symbol = "87-Sr";
	itp->name = "87-Strontium";
	itp->isotope = true;
	itp->M = 86.9088793;
	itp->nsf = 7.41;
	pt_backup.push_back(itp);
	// 88-strontium
	itp = new AtomType(*atp);
	itp->symbol = "88-Sr";
	itp->name = "88-Strontium";
	itp->isotope = true;
	itp->M = 87.9056143;
	itp->nsf = 7.16;
	pt_backup.push_back(itp);
    }
    // yttrium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Y";
	atp->name = "Yttrium";
	atp->z = 39;
	atp->M = 88.905852;
	atp->radius = 1.8237;
	atp->xsf = 39.0;
	atp->nsf = 7.752;
    }
    // zirconium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Zr";
	atp->name = "Zirconium";
	atp->z = 40;
	atp->M = 91.2242;
	atp->radius = 1.6156;
	atp->xsf = 40.0;
	atp->nsf = 7.163;
	// isotopes
	// 90-zirconium
	itp = new AtomType(*atp);
	itp->symbol = "90-Zr";
	itp->name = "90-Zirconium";
	itp->isotope = true;
	itp->M = 89.9047037;
	itp->nsf = 6.5;
	pt_backup.push_back(itp);
	// 91-zirconium
	itp = new AtomType(*atp);
	itp->symbol = "91-Zr";
	itp->name = "91-Zirconium";
	itp->isotope = true;
	itp->M = 90.9056450;
	itp->nsf = 8.8;
	pt_backup.push_back(itp);
	// 92-zirconium
	itp = new AtomType(*atp);
	itp->symbol = "92-Zr";
	itp->name = "92-Zirconium";
	itp->isotope = true;
	itp->M = 91.9050401;
	itp->nsf = 7.5;
	pt_backup.push_back(itp);
	// 94-zirconium
	itp = new AtomType(*atp);
	itp->symbol = "94-Zr";
	itp->name = "94-Zirconium";
	itp->isotope = true;
	itp->M = 93.9063158;
	itp->nsf = 8.3;
	pt_backup.push_back(itp);
	// 96-zirconium
	itp = new AtomType(*atp);
	itp->symbol = "96-Zr";
	itp->name = "96-Zirconium";
	itp->isotope = true;
	itp->M = 95.908276;
	itp->nsf = 5.5;
	pt_backup.push_back(itp);
    }
    // niobium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Nb";
	atp->name = "Niobium";
	atp->z = 41;
	atp->M = 92.906382;
	atp->radius = 1.4318;
	atp->xsf = 41.0;
	atp->nsf = 7.0543;
    }
    // molybdenum
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Mo";
	atp->name = "Molybdenum";
	atp->z = 42;
	atp->M = 95.942;
	atp->radius = 1.3626;
	atp->xsf = 42.0;
	atp->nsf = 6.71520;
	// isotopes
	// 92-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "92-Mo";
	itp->name = "92-Molybdenum";
	itp->isotope = true;
	itp->M = 91.906810;
	itp->nsf = 6.93;
	pt_backup.push_back(itp);
	// 94-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "94-Mo";
	itp->name = "94-Molybdenum";
	itp->isotope = true;
	itp->M = 93.9050876;
	itp->nsf = 6.82;
	pt_backup.push_back(itp);
	// 95-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "95-Mo";
	itp->name = "95-Molybdenum";
	itp->isotope = true;
	itp->M = 94.9058415;
	itp->nsf = 6.93;
	pt_backup.push_back(itp);
	// 96-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "96-Mo";
	itp->name = "96-Molybdenum";
	itp->isotope = true;
	itp->M = 95.9046789;
	itp->nsf = 6.22;
	pt_backup.push_back(itp);
	// 97-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "97-Mo";
	itp->name = "97-Molybdenum";
	itp->isotope = true;
	itp->M = 96.9060210;
	itp->nsf = 7.26;
	pt_backup.push_back(itp);
	// 98-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "98-Mo";
	itp->name = "98-Molybdenum";
	itp->isotope = true;
	itp->M = 97.9054078;
	itp->nsf = 6.60;
	pt_backup.push_back(itp);
	// 100-molybdenum
	itp = new AtomType(*atp);
	itp->symbol = "100-Mo";
	itp->name = "100-Molybdenum";
	itp->isotope = true;
	itp->M = 99.907477;
	itp->nsf = 6.75;
	pt_backup.push_back(itp);
    }
    // technetium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Tc";
	atp->name = "Technetium";
	atp->z = 43;
	atp->M = 98.0;
	atp->radius = 1.3675;
	atp->xsf = 43.0;
	atp->nsf = 6.83;
    }
    // ruthenium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ru";
	atp->name = "Ruthenium";
	atp->z = 44;
	atp->M = 101.072;
	atp->radius = 1.3529;
	atp->xsf = 44.0;
	atp->nsf = 7.022;
    }
    // rhodium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Rh";
	atp->name = "Rhodium";
	atp->z = 45;
	atp->M = 102.905502;
	atp->radius = 1.3450;
	atp->xsf = 45.0;
	atp->nsf = 5.904;
    }
    // palladium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pd";
	atp->name = "Palladium";
	atp->z = 46;
	atp->M = 106.421;
	atp->radius = 1.3755;
	atp->xsf = 46.0;
	atp->nsf = 5.916;
	// isotopes
	// 102-palladium
	itp = new AtomType(*atp);
	itp->symbol = "102-Pd";
	itp->name = "102-Palladium";
	itp->isotope = true;
	itp->M = 101.905608;
	itp->nsf = 7.7;
	pt_backup.push_back(itp);
	// 104-palladium
	itp = new AtomType(*atp);
	itp->symbol = "104-Pd";
	itp->name = "104-Palladium";
	itp->isotope = true;
	itp->M = 103.904035;
	itp->nsf = 7.7;
	pt_backup.push_back(itp);
	// 105-palladium
	itp = new AtomType(*atp);
	itp->symbol = "105-Pd";
	itp->name = "105-Palladium";
	itp->isotope = true;
	itp->M = 104.905084;
	itp->nsf = 5.5;
	pt_backup.push_back(itp);
	// 106-palladium
	itp = new AtomType(*atp);
	itp->symbol = "106-Pd";
	itp->name = "106-Palladium";
	itp->isotope = true;
	itp->M = 105.903483;
	itp->nsf = 6.4;
	pt_backup.push_back(itp);
	// 108-palladium
	itp = new AtomType(*atp);
	itp->symbol = "108-Pd";
	itp->name = "108-Palladium";
	itp->isotope = true;
	itp->M = 107.903894;
	itp->nsf = 4.1;
	pt_backup.push_back(itp);
	// 110-palladium
	itp = new AtomType(*atp);
	itp->symbol = "110-Pd";
	itp->name = "110-Palladium";
	itp->isotope = true;
	itp->M = 109.905152;
	itp->nsf = 7.7;
	pt_backup.push_back(itp);
    }
    // silver
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ag";
	atp->name = "Silver";
	atp->z = 47;
	atp->M = 107.86822;
	atp->radius = 1.4447;
	atp->xsf = 47.0;
	atp->nsf = 5.9227;
	// isotopes
	// 107-silver
	itp = new AtomType(*atp);
	itp->symbol = "107-Ag";
	itp->name = "107-Silver";
	itp->isotope = true;
	itp->M = 106.905093;
	itp->nsf = 7.555;
	pt_backup.push_back(itp);
	// 109-silver
	itp = new AtomType(*atp);
	itp->symbol = "109-Ag";
	itp->name = "109-Silver";
	itp->isotope = true;
	itp->M = 108.904756;
	itp->nsf = 4.165;
	pt_backup.push_back(itp);
    }
    // cadmium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cd";
	atp->name = "Cadmium";
	atp->z = 48;
	atp->M = 112.4118;
	atp->radius = 1.4894;
	atp->xsf = 48.0;
	atp->nsf = 4.835;
	// isotopes
	// 106-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "106-Cd";
	itp->name = "106-Cadmium";
	itp->isotope = true;
	itp->M = 105.906458;
	itp->nsf = 5.0;
	pt_backup.push_back(itp);
	// 108-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "108-Cd";
	itp->name = "108-Cadmium";
	itp->isotope = true;
	itp->M = 107.904183;
	itp->nsf = 5.31;
	pt_backup.push_back(itp);
	// 110-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "110-Cd";
	itp->name = "110-Cadmium";
	itp->isotope = true;
	itp->M = 109.903006;
	itp->nsf = 5.78;
	pt_backup.push_back(itp);
	// 111-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "111-Cd";
	itp->name = "111-Cadmium";
	itp->isotope = true;
	itp->M = 110.904182;
	itp->nsf = 6.47;
	pt_backup.push_back(itp);
	// 112-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "112-Cd";
	itp->name = "112-Cadmium";
	itp->isotope = true;
	itp->M = 111.9027572;
	itp->nsf = 6.34;
	pt_backup.push_back(itp);
	// 113-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "113-Cd";
	itp->name = "113-Cadmium";
	itp->isotope = true;
	itp->M = 112.9044009;
	itp->nsf = -8.0;
	pt_backup.push_back(itp);
	// 114-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "114-Cd";
	itp->name = "114-Cadmium";
	itp->isotope = true;
	itp->M = 113.9033581;
	itp->nsf = 7.48;
	pt_backup.push_back(itp);
	// 116-cadmium
	itp = new AtomType(*atp);
	itp->symbol = "116-Cd";
	itp->name = "116-Cadmium";
	itp->isotope = true;
	itp->M = 115.904755;
	itp->nsf = 6.26;
	pt_backup.push_back(itp);
    }
    // indium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "In";
	atp->name = "Indium";
	atp->z = 49;
	atp->M = 114.8183;
	atp->radius = 1.6662;
	atp->xsf = 49.0;
	atp->nsf = 4.0652;
	// isotopes
	// 113-indium
	itp = new AtomType(*atp);
	itp->symbol = "113-In";
	itp->name = "113-Indium";
	itp->isotope = true;
	itp->M = 112.904061;
	itp->nsf = 5.39;
	pt_backup.push_back(itp);
	// 115-indium
	itp = new AtomType(*atp);
	itp->symbol = "115-In";
	itp->name = "115-Indium";
	itp->isotope = true;
	itp->M = 114.903878;
	itp->nsf = 4.00;
	pt_backup.push_back(itp);
    }
    // tin
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sn";
	atp->name = "Tin";
	atp->z = 50;
	atp->M = 118.7107;
	atp->radius = 1.5375;
	atp->xsf = 50.0;
	atp->nsf = 6.2252;
	// isotopes
	// 112-tin
	itp = new AtomType(*atp);
	itp->symbol = "112-Sn";
	itp->name = "112-Tin";
	itp->isotope = true;
	itp->M = 111.904821;
	itp->nsf = 6.0;
	pt_backup.push_back(itp);
	// 114-tin
	itp = new AtomType(*atp);
	itp->symbol = "114-Sn";
	itp->name = "114-Tin";
	itp->isotope = true;
	itp->M = 113.902782;
	itp->nsf = 6.0;
	pt_backup.push_back(itp);
	// 115-tin
	itp = new AtomType(*atp);
	itp->symbol = "115-Sn";
	itp->name = "115-Tin";
	itp->isotope = true;
	itp->M = 114.903346;
	itp->nsf = 6.0;
	pt_backup.push_back(itp);
	// 116-tin
	itp = new AtomType(*atp);
	itp->symbol = "116-Sn";
	itp->name = "116-Tin";
	itp->isotope = true;
	itp->M = 115.901744;
	itp->nsf = 6.1;
	pt_backup.push_back(itp);
	// 117-tin
	itp = new AtomType(*atp);
	itp->symbol = "117-Sn";
	itp->name = "117-Tin";
	itp->isotope = true;
	itp->M = 116.902954;
	itp->nsf = 6.59;
	pt_backup.push_back(itp);
	// 118-tin
	itp = new AtomType(*atp);
	itp->symbol = "118-Sn";
	itp->name = "118-Tin";
	itp->isotope = true;
	itp->M = 117.901606;
	itp->nsf = 6.23;
	pt_backup.push_back(itp);
	// 119-tin
	itp = new AtomType(*atp);
	itp->symbol = "119-Sn";
	itp->name = "119-Tin";
	itp->isotope = true;
	itp->M = 118.903309;
	itp->nsf = 6.28;
	pt_backup.push_back(itp);
	// 120-tin
	itp = new AtomType(*atp);
	itp->symbol = "120-Sn";
	itp->name = "120-Tin";
	itp->isotope = true;
	itp->M = 119.9021966;
	itp->nsf = 6.67;
	pt_backup.push_back(itp);
	// 122-tin
	itp = new AtomType(*atp);
	itp->symbol = "122-Sn";
	itp->name = "122-Tin";
	itp->isotope = true;
	itp->M = 121.9034401;
	itp->nsf = 5.93;
	pt_backup.push_back(itp);
	// 124-tin
	itp = new AtomType(*atp);
	itp->symbol = "124-Sn";
	itp->name = "124-Tin";
	itp->isotope = true;
	itp->M = 123.9052746;
	itp->nsf = 6.15;
	pt_backup.push_back(itp);
    }
    // antimony
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sb";
	atp->name = "Antimony";
	atp->z = 51;
	atp->M = 121.7601;
	atp->radius = 1.4000;
	atp->xsf = 51.0;
	atp->nsf = 5.573;
	// isotopes
	// 121-antimony
	itp = new AtomType(*atp);
	itp->symbol = "121-Sb";
	itp->name = "121-Antimony";
	itp->isotope = true;
	itp->M = 120.9038180;
	itp->nsf = 5.71;
	pt_backup.push_back(itp);
	// 123-antimony
	itp = new AtomType(*atp);
	itp->symbol = "123-Sb";
	itp->name = "123-Antimony";
	itp->isotope = true;
	itp->M = 122.9042157;
	itp->nsf = 5.38;
	pt_backup.push_back(itp);
    }
    // tellurium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Te";
	atp->name = "Tellurium";
	atp->z = 52;
	atp->M = 127.603;
	atp->radius = 1.3600;
	atp->xsf = 52.0;
	atp->nsf = 5.682;
	// isotopes
	// 120-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "120-Te";
	itp->name = "120-Tellurium";
	itp->isotope = true;
	itp->M = 119.904020;
	itp->nsf = 5.3;
	pt_backup.push_back(itp);
	// 122-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "122-Te";
	itp->name = "122-Tellurium";
	itp->isotope = true;
	itp->M = 121.9030471;
	itp->nsf = 3.8;
	pt_backup.push_back(itp);
	// 123-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "123-Te";
	itp->name = "123-Tellurium";
	itp->isotope = true;
	itp->M = 122.9042730;
	itp->nsf = -0.05;
	pt_backup.push_back(itp);
	// 124-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "124-Te";
	itp->name = "124-Tellurium";
	itp->isotope = true;
	itp->M = 123.9028195;
	itp->nsf = 7.95;
	pt_backup.push_back(itp);
	// 125-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "125-Te";
	itp->name = "125-Tellurium";
	itp->isotope = true;
	itp->M = 124.9044247;
	itp->nsf = 5.01;
	pt_backup.push_back(itp);
	// 126-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "126-Te";
	itp->name = "126-Tellurium";
	itp->isotope = true;
	itp->M = 125.9033055;
	itp->nsf = 5.55;
	pt_backup.push_back(itp);
	// 128-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "128-Te";
	itp->name = "128-Tellurium";
	itp->isotope = true;
	itp->M = 127.9044614;
	itp->nsf = 5.88;
	pt_backup.push_back(itp);
	// 130-tellurium
	itp = new AtomType(*atp);
	itp->symbol = "130-Te";
	itp->name = "130-Tellurium";
	itp->isotope = true;
	itp->M = 129.9062228;
	itp->nsf = 6.01;
	pt_backup.push_back(itp);
    }
    // iodine
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "I";
	atp->name = "Iodine";
	atp->z = 53;
	atp->M = 126.904473;
	atp->radius = 1.3300;
	atp->xsf = 53.0;
	atp->nsf = 5.282;
    }
    // xenon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Xe";
	atp->name = "Xenon";
	atp->z = 54;
	atp->M = 131.2936;
	atp->radius = 2.2000;
	atp->xsf = 54.0;
	atp->nsf = 4.694;
    }
    // cesium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cs";
	atp->name = "Cesium";
	atp->z = 55;
	atp->M = 132.905452;
	atp->radius = 2.6325;
	atp->xsf = 55.0;
	atp->nsf = 5.422;
    }
    // barium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ba";
	atp->name = "Barium";
	atp->z = 56;
	atp->M = 137.3277;
	atp->radius = 2.1705;
	atp->xsf = 56.0;
	atp->nsf = 5.073;
	// isotopes
	// 130-barium
	itp = new AtomType(*atp);
	itp->symbol = "130-Ba";
	itp->name = "130-Barium";
	itp->isotope = true;
	itp->M = 129.906310;
	itp->nsf = -3.6;
	pt_backup.push_back(itp);
	// 132-barium
	itp = new AtomType(*atp);
	itp->symbol = "132-Ba";
	itp->name = "132-Barium";
	itp->isotope = true;
	itp->M = 131.905056;
	itp->nsf = 7.8;
	pt_backup.push_back(itp);
	// 134-barium
	itp = new AtomType(*atp);
	itp->symbol = "134-Ba";
	itp->name = "134-Barium";
	itp->isotope = true;
	itp->M = 133.904503;
	itp->nsf = 5.7;
	pt_backup.push_back(itp);
	// 135-barium
	itp = new AtomType(*atp);
	itp->symbol = "135-Ba";
	itp->name = "135-Barium";
	itp->isotope = true;
	itp->M = 134.905683;
	itp->nsf = 4.66;
	pt_backup.push_back(itp);
	// 136-barium
	itp = new AtomType(*atp);
	itp->symbol = "136-Ba";
	itp->name = "136-Barium";
	itp->isotope = true;
	itp->M = 135.904570;
	itp->nsf = 4.90;
	pt_backup.push_back(itp);
	// 137-barium
	itp = new AtomType(*atp);
	itp->symbol = "137-Ba";
	itp->name = "137-Barium";
	itp->isotope = true;
	itp->M = 136.905821;
	itp->nsf = 6.82;
	pt_backup.push_back(itp);
	// 138-barium
	itp = new AtomType(*atp);
	itp->symbol = "138-Ba";
	itp->name = "138-Barium";
	itp->isotope = true;
	itp->M = 137.905241;
	itp->nsf = 4.83;
	pt_backup.push_back(itp);
    }
    // lanthanum
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "La";
	atp->name = "Lanthanum";
	atp->z = 57;
	atp->M = 138.90552;
	atp->radius = 1.8725;
	atp->xsf = 57.0;
	atp->nsf = 8.244;
	// isotopes
	// 138-lanthanum
	itp = new AtomType(*atp);
	itp->symbol = "138-La";
	itp->name = "138-Lanthanum";
	itp->isotope = true;
	itp->M = 137.907107;
	itp->nsf = 8.0;
	pt_backup.push_back(itp);
	// 139-lanthanum
	itp = new AtomType(*atp);
	itp->symbol = "139-La";
	itp->name = "139-Lanthanum";
	itp->isotope = true;
	itp->M = 138.906348;
	itp->nsf = 8.24;
	pt_backup.push_back(itp);
    }
    // cerium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ce";
	atp->name = "Cerium";
	atp->z = 58;
	atp->M = 140.1161;
	atp->radius = 1.8243;
	atp->xsf = 58.0;
	atp->nsf = 4.842;
	// isotopes
	// 136-cerium
	itp = new AtomType(*atp);
	itp->symbol = "136-Ce";
	itp->name = "136-Cerium";
	itp->isotope = true;
	itp->M = 135.907140;
	itp->nsf = 5.76;
	pt_backup.push_back(itp);
	// 138-cerium
	itp = new AtomType(*atp);
	itp->symbol = "138-Ce";
	itp->name = "138-Cerium";
	itp->isotope = true;
	itp->M = 137.905986;
	itp->nsf = 6.65;
	pt_backup.push_back(itp);
	// 140-cerium
	itp = new AtomType(*atp);
	itp->symbol = "140-Ce";
	itp->name = "140-Cerium";
	itp->isotope = true;
	itp->M = 139.905434;
	itp->nsf = 4.81;
	pt_backup.push_back(itp);
	// 142-cerium
	itp = new AtomType(*atp);
	itp->symbol = "142-Ce";
	itp->name = "142-Cerium";
	itp->isotope = true;
	itp->M = 141.909240;
	itp->nsf = 4.72;
	pt_backup.push_back(itp);
    }
    // praseodymium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pr";
	atp->name = "Praseodymium";
	atp->z = 59;
	atp->M = 140.907652;
	atp->radius = 1.8362;
	atp->xsf = 59.0;
	atp->nsf = 4.585;
    }
    // neodymium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Nd";
	atp->name = "Neodymium";
	atp->z = 60;
	atp->M = 144.243;
	atp->radius = 1.8295;
	atp->xsf = 60.0;
	atp->nsf = 7.695;
	// isotopes
	// 142-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "142-Nd";
	itp->name = "142-Neodymium";
	itp->isotope = true;
	itp->M = 141.907719;
	itp->nsf = 7.7;
	pt_backup.push_back(itp);
	// 143-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "143-Nd";
	itp->name = "143-Neodymium";
	itp->isotope = true;
	itp->M = 142.909810;
	itp->nsf = 14.0;
	pt_backup.push_back(itp);
	// 144-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "144-Nd";
	itp->name = "144-Neodymium";
	itp->isotope = true;
	itp->M = 143.910083;
	itp->nsf = 2.8;
	pt_backup.push_back(itp);
	// 145-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "145-Nd";
	itp->name = "145-Neodymium";
	itp->isotope = true;
	itp->M = 144.912569;
	itp->nsf = 14.0;
	pt_backup.push_back(itp);
	// 146-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "146-Nd";
	itp->name = "146-Neodymium";
	itp->isotope = true;
	itp->M = 145.913112;
	itp->nsf = 8.7;
	pt_backup.push_back(itp);
	// 148-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "148-Nd";
	itp->name = "148-Neodymium";
	itp->isotope = true;
	itp->M = 147.916889;
	itp->nsf = 5.7;
	pt_backup.push_back(itp);
	// 150-neodymium
	itp = new AtomType(*atp);
	itp->symbol = "150-Nd";
	itp->name = "150-Neodymium";
	itp->isotope = true;
	itp->M = 149.920887;
	itp->nsf = 5.28;
	pt_backup.push_back(itp);
    }
    // promethium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pm";
	atp->name = "Promethium";
	atp->z = 61;
	atp->M = 145.0;
	atp->radius = 1.8090;
	atp->xsf = 61.0;
	atp->nsf = 12.64;
    }
    // samarium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sm";
	atp->name = "Samarium";
	atp->z = 62;
	atp->M = 150.363;
	atp->radius = 1.8040;
	atp->xsf = 62.0;
	atp->nsf = 0.005;
	// isotopes
	// 144-samarium
	itp = new AtomType(*atp);
	itp->symbol = "144-Sm";
	itp->name = "144-Samarium";
	itp->isotope = true;
	itp->M = 143.911995;
	itp->nsf = -3.0;
	pt_backup.push_back(itp);
	// 147-samarium
	itp = new AtomType(*atp);
	itp->symbol = "147-Sm";
	itp->name = "147-Samarium";
	itp->isotope = true;
	itp->M = 146.914893;
	itp->nsf = 14.0;
	pt_backup.push_back(itp);
	// 148-samarium
	itp = new AtomType(*atp);
	itp->symbol = "148-Sm";
	itp->name = "148-Samarium";
	itp->isotope = true;
	itp->M = 147.914818;
	itp->nsf = -3.0;
	pt_backup.push_back(itp);
	// 149-samarium
	itp = new AtomType(*atp);
	itp->symbol = "149-Sm";
	itp->name = "149-Samarium";
	itp->isotope = true;
	itp->M = 148.917180;
	itp->nsf = 18.7;
	pt_backup.push_back(itp);
	// 150-samarium
	itp = new AtomType(*atp);
	itp->symbol = "150-Sm";
	itp->name = "150-Samarium";
	itp->isotope = true;
	itp->M = 149.917271;
	itp->nsf = 14.0;
	pt_backup.push_back(itp);
	// 152-samarium
	itp = new AtomType(*atp);
	itp->symbol = "152-Sm";
	itp->name = "152-Samarium";
	itp->isotope = true;
	itp->M = 151.919728;
	itp->nsf = -5.0;
	pt_backup.push_back(itp);
	// 154-samarium
	itp = new AtomType(*atp);
	itp->symbol = "154-Sm";
	itp->name = "154-Samarium";
	itp->isotope = true;
	itp->M = 153.922205;
	itp->nsf = 8.0;
	pt_backup.push_back(itp);
    }
    // europium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Eu";
	atp->name = "Europium";
	atp->z = 63;
	atp->M = 151.9641;
	atp->radius = 1.9840;
	atp->xsf = 63.0;
	atp->nsf = 5.33;
	// isotopes
	// 153-europium
	itp = new AtomType(*atp);
	itp->symbol = "153-Eu";
	itp->name = "153-Europium";
	itp->isotope = true;
	itp->M = 152.921226;
	itp->nsf = 8.22;
	pt_backup.push_back(itp);
    }
    // gadolinium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Gd";
	atp->name = "Gadolinium";
	atp->z = 64;
	atp->M = 157.253;
	atp->radius = 1.8180;
	atp->xsf = 64.0;
	atp->nsf = 9.52;
	// isotopes
	// 152-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "152-Gd";
	itp->name = "152-Gadolinium";
	itp->isotope = true;
	itp->M = 151.919788;
	itp->nsf = 10.0;
	pt_backup.push_back(itp);
	// 154-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "154-Gd";
	itp->name = "154-Gadolinium";
	itp->isotope = true;
	itp->M = 153.920862;
	itp->nsf = 10.0;
	pt_backup.push_back(itp);
	// 155-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "155-Gd";
	itp->name = "155-Gadolinium";
	itp->isotope = true;
	itp->M = 154.922619;
	itp->nsf = 13.8;
	pt_backup.push_back(itp);
	// 156-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "156-Gd";
	itp->name = "156-Gadolinium";
	itp->isotope = true;
	itp->M = 155.922120;
	itp->nsf = 6.3;
	pt_backup.push_back(itp);
	// 157-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "157-Gd";
	itp->name = "157-Gadolinium";
	itp->isotope = true;
	itp->M = 156.923957;
	itp->nsf = 4.0;
	pt_backup.push_back(itp);
	// 158-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "158-Gd";
	itp->name = "158-Gadolinium";
	itp->isotope = true;
	itp->M = 157.924101;
	itp->nsf = 9.0;
	pt_backup.push_back(itp);
	// 160-gadolinium
	itp = new AtomType(*atp);
	itp->symbol = "160-Gd";
	itp->name = "160-Gadolinium";
	itp->isotope = true;
	itp->M = 159.927051;
	itp->nsf = 9.15;
	pt_backup.push_back(itp);
    }
    // terbium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Tb";
	atp->name = "Terbium";
	atp->z = 65;
	atp->M = 158.925342;
	atp->radius = 1.8005;
	atp->xsf = 65.0;
	atp->nsf = 7.342;
    }
    // dysprosium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Dy";
	atp->name = "Dysprosium";
	atp->z = 66;
	atp->M = 162.5001;
	atp->radius = 1.7951;
	atp->xsf = 66.0;
	atp->nsf = 16.93;
	// isotopes
	// 156-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "156-Dy";
	itp->name = "156-Dysprosium";
	itp->isotope = true;
	itp->M = 155.924278;
	itp->nsf = 6.1;
	pt_backup.push_back(itp);
	// 158-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "158-Dy";
	itp->name = "158-Dysprosium";
	itp->isotope = true;
	itp->M = 157.924405;
	itp->nsf = 6.0;
	pt_backup.push_back(itp);
	// 160-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "160-Dy";
	itp->name = "160-Dysprosium";
	itp->isotope = true;
	itp->M = 159.925194;
	itp->nsf = 6.7;
	pt_backup.push_back(itp);
	// 161-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "161-Dy";
	itp->name = "161-Dysprosium";
	itp->isotope = true;
	itp->M = 160.926930;
	itp->nsf = 10.3;
	pt_backup.push_back(itp);
	// 162-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "162-Dy";
	itp->name = "162-Dysprosium";
	itp->isotope = true;
	itp->M = 161.926795;
	itp->nsf = -1.4;
	pt_backup.push_back(itp);
	// 163-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "163-Dy";
	itp->name = "163-Dysprosium";
	itp->isotope = true;
	itp->M = 162.928728;
	itp->nsf = 5.0;
	pt_backup.push_back(itp);
	// 164-dysprosium
	itp = new AtomType(*atp);
	itp->symbol = "164-Dy";
	itp->name = "164-Dysprosium";
	itp->isotope = true;
	itp->M = 163.929171;
	itp->nsf = 49.4;
	pt_backup.push_back(itp);
    }
    // holmium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ho";
	atp->name = "Holmium";
	atp->z = 67;
	atp->M = 164.930322;
	atp->radius = 1.7886;
	atp->xsf = 67.0;
	atp->nsf = 8.443;
    }
    // erbium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Er";
	atp->name = "Erbium";
	atp->z = 68;
	atp->M = 167.2593;
	atp->radius = 1.7794;
	atp->xsf = 68.0;
	atp->nsf = 7.792;
	// isotopes
	// 162-erbium
	itp = new AtomType(*atp);
	itp->symbol = "162-Er";
	itp->name = "162-Erbium";
	itp->isotope = true;
	itp->M = 161.928775;
	itp->nsf = 9.01;
	pt_backup.push_back(itp);
	// 164-erbium
	itp = new AtomType(*atp);
	itp->symbol = "164-Er";
	itp->name = "164-Erbium";
	itp->isotope = true;
	itp->M = 163.929197;
	itp->nsf = 7.95;
	pt_backup.push_back(itp);
	// 166-erbium
	itp = new AtomType(*atp);
	itp->symbol = "166-Er";
	itp->name = "166-Erbium";
	itp->isotope = true;
	itp->M = 165.930290;
	itp->nsf = 10.51;
	pt_backup.push_back(itp);
	// 167-erbium
	itp = new AtomType(*atp);
	itp->symbol = "167-Er";
	itp->name = "167-Erbium";
	itp->isotope = true;
	itp->M = 166.932045;
	itp->nsf = 3.06;
	pt_backup.push_back(itp);
	// 168-erbium
	itp = new AtomType(*atp);
	itp->symbol = "168-Er";
	itp->name = "168-Erbium";
	itp->isotope = true;
	itp->M = 167.932368;
	itp->nsf = 7.43;
	pt_backup.push_back(itp);
	// 170-erbium
	itp = new AtomType(*atp);
	itp->symbol = "170-Er";
	itp->name = "170-Erbium";
	itp->isotope = true;
	itp->M = 169.935460;
	itp->nsf = 9.61;
	pt_backup.push_back(itp);
    }
    // thulium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Tm";
	atp->name = "Thulium";
	atp->z = 69;
	atp->M = 168.934212;
	atp->radius = 1.7687;
	atp->xsf = 69.0;
	atp->nsf = 7.073;
    }
    // ytterbium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Yb";
	atp->name = "Ytterbium";
	atp->z = 70;
	atp->M = 173.043;
	atp->radius = 1.9396;
	atp->xsf = 70.0;
	atp->nsf = 12.413;
	// isotopes
	// 168-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "168-Yb";
	itp->name = "168-Ytterbium";
	itp->isotope = true;
	itp->M = 167.933894;
	itp->nsf = -4.07;
	pt_backup.push_back(itp);
	// 170-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "170-Yb";
	itp->name = "170-Ytterbium";
	itp->isotope = true;
	itp->M = 169.934759;
	itp->nsf = 6.8;
	pt_backup.push_back(itp);
	// 171-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "171-Yb";
	itp->name = "171-Ytterbium";
	itp->isotope = true;
	itp->M = 170.936322;
	itp->nsf = 9.7;
	pt_backup.push_back(itp);
	// 172-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "172-Yb";
	itp->name = "172-Ytterbium";
	itp->isotope = true;
	itp->M = 171.9363777;
	itp->nsf = 9.5;
	pt_backup.push_back(itp);
	// 173-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "173-Yb";
	itp->name = "173-Ytterbium";
	itp->isotope = true;
	itp->M = 172.9382068;
	itp->nsf = 9.56;
	pt_backup.push_back(itp);
	// 174-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "174-Yb";
	itp->name = "174-Ytterbium";
	itp->isotope = true;
	itp->M = 173.9388581;
	itp->nsf = 19.2;
	pt_backup.push_back(itp);
	// 176-ytterbium
	itp = new AtomType(*atp);
	itp->symbol = "176-Yb";
	itp->name = "176-Ytterbium";
	itp->isotope = true;
	itp->M = 175.942568;
	itp->nsf = 8.7;
	pt_backup.push_back(itp);
    }
    // lutetium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Lu";
	atp->name = "Lutetium";
	atp->z = 71;
	atp->M = 174.9671;
	atp->radius = 1.7515;
	atp->xsf = 71.0;
	atp->nsf = 7.213;
	// isotopes
	// 175-lutetium
	itp = new AtomType(*atp);
	itp->symbol = "175-Lu";
	itp->name = "175-Lutetium";
	itp->isotope = true;
	itp->M = 174.9407679;
	itp->nsf = 7.28;
	pt_backup.push_back(itp);
	// 176-lutetium
	itp = new AtomType(*atp);
	itp->symbol = "176-Lu";
	itp->name = "176-Lutetium";
	itp->isotope = true;
	itp->M = 175.9426824;
	itp->nsf = 6.1;
	pt_backup.push_back(itp);
    }
    // hafnium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Hf";
	atp->name = "Hafnium";
	atp->z = 72;
	atp->M = 178.492;
	atp->radius = 1.5973;
	atp->xsf = 72.0;
	atp->nsf = 7.7714;
	// isotopes
	// 174-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "174-Hf";
	itp->name = "174-Hafnium";
	itp->isotope = true;
	itp->M = 173.940040;
	itp->nsf = 10.9;
	pt_backup.push_back(itp);
	// 176-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "176-Hf";
	itp->name = "176-Hafnium";
	itp->isotope = true;
	itp->M = 175.9414018;
	itp->nsf = 6.61;
	pt_backup.push_back(itp);
	// 177-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "177-Hf";
	itp->name = "177-Hafnium";
	itp->isotope = true;
	itp->M = 176.9432200;
	itp->nsf = 0.8;
	pt_backup.push_back(itp);
	// 178-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "178-Hf";
	itp->name = "178-Hafnium";
	itp->isotope = true;
	itp->M = 177.9436977;
	itp->nsf = 5.9;
	pt_backup.push_back(itp);
	// 179-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "179-Hf";
	itp->name = "179-Hafnium";
	itp->isotope = true;
	itp->M = 178.9458151;
	itp->nsf = 7.46;
	pt_backup.push_back(itp);
	// 180-hafnium
	itp = new AtomType(*atp);
	itp->symbol = "180-Hf";
	itp->name = "180-Hafnium";
	itp->isotope = true;
	itp->M = 179.9465488;
	itp->nsf = 13.2;
	pt_backup.push_back(itp);
    }
    // tantalum
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ta";
	atp->name = "Tantalum";
	atp->z = 73;
	atp->M = 180.94791;
	atp->radius = 1.4280;
	atp->xsf = 73.0;
	atp->nsf = 6.917;
	// isotopes
	// 180-tantalum
	itp = new AtomType(*atp);
	itp->symbol = "180-Ta";
	itp->name = "180-Tantalum";
	itp->isotope = true;
	itp->M = 179.947466;
	itp->nsf = 7.0;
	pt_backup.push_back(itp);
	// 181-tantalum
	itp = new AtomType(*atp);
	itp->symbol = "181-Ta";
	itp->name = "181-Tantalum";
	itp->isotope = true;
	itp->M = 180.947996;
	itp->nsf = 6.91;
	pt_backup.push_back(itp);
    }
    // tungsten
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "W";
	atp->name = "Tungsten";
	atp->z = 74;
	atp->M = 183.841;
	atp->radius = 1.3705;
	atp->xsf = 74.0;
	atp->nsf = 4.75518;
	// isotopes
	// 180-tungsten
	itp = new AtomType(*atp);
	itp->symbol = "180-W";
	itp->name = "180-Tungsten";
	itp->isotope = true;
	itp->M = 179.946706;
	itp->nsf = 5.0;
	pt_backup.push_back(itp);
	// 182-tungsten
	itp = new AtomType(*atp);
	itp->symbol = "182-W";
	itp->name = "182-Tungsten";
	itp->isotope = true;
	itp->M = 181.948206;
	itp->nsf = 7.04;
	pt_backup.push_back(itp);
	// 183-tungsten
	itp = new AtomType(*atp);
	itp->symbol = "183-W";
	itp->name = "183-Tungsten";
	itp->isotope = true;
	itp->M = 182.9502245;
	itp->nsf = 6.59;
	pt_backup.push_back(itp);
	// 184-tungsten
	itp = new AtomType(*atp);
	itp->symbol = "184-W";
	itp->name = "184-Tungsten";
	itp->isotope = true;
	itp->M = 183.9509326;
	itp->nsf = 7.55;
	pt_backup.push_back(itp);
	// 186-tungsten
	itp = new AtomType(*atp);
	itp->symbol = "186-W";
	itp->name = "186-Tungsten";
	itp->isotope = true;
	itp->M = 185.954362;
	itp->nsf = -0.73;
	pt_backup.push_back(itp);
    }
    // rhenium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Re";
	atp->name = "Rhenium";
	atp->z = 75;
	atp->M = 186.2071;
	atp->radius = 1.3800;
	atp->xsf = 75.0;
	atp->nsf = 9.22;
	// isotopes
	// 185-rhenium
	itp = new AtomType(*atp);
	itp->symbol = "185-Re";
	itp->name = "185-Rhenium";
	itp->isotope = true;
	itp->M = 184.9529557;
	itp->nsf = 9.0;
	pt_backup.push_back(itp);
	// 187-rhenium
	itp = new AtomType(*atp);
	itp->symbol = "187-Re";
	itp->name = "187-Rhenium";
	itp->isotope = true;
	itp->M = 186.9557508;
	itp->nsf = 9.3;
	pt_backup.push_back(itp);
    }
    // osmium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Os";
	atp->name = "Osmium";
	atp->z = 76;
	atp->M = 190.233;
	atp->radius = 1.3676;
	atp->xsf = 76.0;
	atp->nsf = 10.72;
	// isotopes
	// 184-osmium
	itp = new AtomType(*atp);
	itp->symbol = "184-Os";
	itp->name = "184-Osmium";
	itp->isotope = true;
	itp->M = 183.952491;
	itp->nsf = 10.0;
	pt_backup.push_back(itp);
	// 186-osmium
	itp = new AtomType(*atp);
	itp->symbol = "186-Os";
	itp->name = "186-Osmium";
	itp->isotope = true;
	itp->M = 185.953838;
	itp->nsf = 12.0;
	pt_backup.push_back(itp);
	// 187-osmium
	itp = new AtomType(*atp);
	itp->symbol = "187-Os";
	itp->name = "187-Osmium";
	itp->isotope = true;
	itp->M = 186.9557479;
	itp->nsf = 10.0;
	pt_backup.push_back(itp);
	// 188-osmium
	itp = new AtomType(*atp);
	itp->symbol = "188-Os";
	itp->name = "188-Osmium";
	itp->isotope = true;
	itp->M = 187.9558360;
	itp->nsf = 7.8;
	pt_backup.push_back(itp);
	// 189-osmium
	itp = new AtomType(*atp);
	itp->symbol = "189-Os";
	itp->name = "189-Osmium";
	itp->isotope = true;
	itp->M = 188.9581449;
	itp->nsf = 11.0;
	pt_backup.push_back(itp);
	// 190-osmium
	itp = new AtomType(*atp);
	itp->symbol = "190-Os";
	itp->name = "190-Osmium";
	itp->isotope = true;
	itp->M = 189.958445;
	itp->nsf = 11.4;
	pt_backup.push_back(itp);
	// 192-osmium
	itp = new AtomType(*atp);
	itp->symbol = "192-Os";
	itp->name = "192-Osmium";
	itp->isotope = true;
	itp->M = 191.961479;
	itp->nsf = 11.9;
	pt_backup.push_back(itp);
    }
    // iridium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ir";
	atp->name = "Iridium";
	atp->z = 77;
	atp->M = 192.2173;
	atp->radius = 1.3573;
	atp->xsf = 77.0;
	atp->nsf = 10.63;
    }
    // platinum
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pt";
	atp->name = "Platinum";
	atp->z = 78;
	atp->M = 195.0782;
	atp->radius = 1.3873;
	atp->xsf = 78.0;
	atp->nsf = 9.601;
	// isotopes
	// 190-platinum
	itp = new AtomType(*atp);
	itp->symbol = "190-Pt";
	itp->name = "190-Platinum";
	itp->isotope = true;
	itp->M = 189.959930;
	itp->nsf = 9.0;
	pt_backup.push_back(itp);
	// 192-platinum
	itp = new AtomType(*atp);
	itp->symbol = "192-Pt";
	itp->name = "192-Platinum";
	itp->isotope = true;
	itp->M = 191.961035;
	itp->nsf = 9.9;
	pt_backup.push_back(itp);
	// 194-platinum
	itp = new AtomType(*atp);
	itp->symbol = "194-Pt";
	itp->name = "194-Platinum";
	itp->isotope = true;
	itp->M = 193.962664;
	itp->nsf = 10.55;
	pt_backup.push_back(itp);
	// 195-platinum
	itp = new AtomType(*atp);
	itp->symbol = "195-Pt";
	itp->name = "195-Platinum";
	itp->isotope = true;
	itp->M = 194.964774;
	itp->nsf = 8.91;
	pt_backup.push_back(itp);
	// 196-platinum
	itp = new AtomType(*atp);
	itp->symbol = "196-Pt";
	itp->name = "196-Platinum";
	itp->isotope = true;
	itp->M = 195.964935;
	itp->nsf = 9.89;
	pt_backup.push_back(itp);
	// 198-platinum
	itp = new AtomType(*atp);
	itp->symbol = "198-Pt";
	itp->name = "198-Platinum";
	itp->isotope = true;
	itp->M = 197.967876;
	itp->nsf = 7.8;
	pt_backup.push_back(itp);
    }
    // gold
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Au";
	atp->name = "Gold";
	atp->z = 79;
	atp->M = 196.966552;
	atp->radius = 1.4419;
	atp->xsf = 79.0;
	atp->nsf = 7.907;
    }
    // mercury
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Hg";
	atp->name = "Mercury";
	atp->z = 80;
	atp->M = 200.592;
	atp->radius = 1.5025;
	atp->xsf = 80.0;
	atp->nsf = 12.59545;
	// isotopes
	// 196-mercury
	itp = new AtomType(*atp);
	itp->symbol = "196-Hg";
	itp->name = "196-Mercury";
	itp->isotope = true;
	itp->M = 195.965815;
	itp->nsf = 30.3;
	pt_backup.push_back(itp);
	// 199-mercury
	itp = new AtomType(*atp);
	itp->symbol = "199-Hg";
	itp->name = "199-Mercury";
	itp->isotope = true;
	itp->M = 198.968262;
	itp->nsf = 16.9;
	pt_backup.push_back(itp);
	// 202-mercury
	itp = new AtomType(*atp);
	itp->symbol = "202-Hg";
	itp->name = "202-Mercury";
	itp->isotope = true;
	itp->M = 201.970626;
	itp->nsf = 11.002;
	pt_backup.push_back(itp);
    }
    // thallium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Tl";
	atp->name = "Thallium";
	atp->z = 81;
	atp->M = 204.38332;
	atp->radius = 1.7283;
	atp->xsf = 81.0;
	atp->nsf = 8.7765;
	// isotopes
	// 203-thalium
	itp = new AtomType(*atp);
	itp->symbol = "203-Tl";
	itp->name = "203-Thallium";
	itp->isotope = true;
	itp->M = 202.972329;
	itp->nsf = 8.51;
	pt_backup.push_back(itp);
	// 205-thalium
	itp = new AtomType(*atp);
	itp->symbol = "205-Tl";
	itp->name = "205-Thallium";
	itp->isotope = true;
	itp->M = 204.974412;
	itp->nsf = 8.87;
	pt_backup.push_back(itp);
    }
    // lead
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pb";
	atp->name = "Lead";
	atp->z = 82;
	atp->M = 207.21;
	atp->radius = 1.7501;
	atp->xsf = 82.0;
	atp->nsf = 9.4012;
	// isotopes
	// 204-lead
	itp = new AtomType(*atp);
	itp->symbol = "204-Pb";
	itp->name = "204-Lead";
	itp->isotope = true;
	itp->M = 203.973029;
	itp->nsf = 10.893;
	pt_backup.push_back(itp);
	// 206-lead
	itp = new AtomType(*atp);
	itp->symbol = "206-Pb";
	itp->name = "206-Lead";
	itp->isotope = true;
	itp->M = 205.974449;
	itp->nsf = 9.2221;
	pt_backup.push_back(itp);
	// 207-lead
	itp = new AtomType(*atp);
	itp->symbol = "207-Pb";
	itp->name = "207-Lead";
	itp->isotope = true;
	itp->M = 206.975881;
	itp->nsf = 9.286;
	pt_backup.push_back(itp);
	// 208-lead
	itp = new AtomType(*atp);
	itp->symbol = "208-Pb";
	itp->name = "208-Lead";
	itp->isotope = true;
	itp->M = 207.976636;
	itp->nsf = 9.494;
	pt_backup.push_back(itp);
    }
    // bismuth
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Bi";
	atp->name = "Bismuth";
	atp->z = 83;
	atp->M = 208.980382;
	atp->radius = 1.4600;
	atp->xsf = 83.0;
	atp->nsf = 8.5322;
    }
    // polonium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Po";
	atp->name = "Polonium";
	atp->z = 84;
	atp->M = 209.0;
	atp->radius = 1.4600;
	atp->xsf = 84.0;
	atp->nsf = 0.0;
    }
    // astatine
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "At";
	atp->name = "Astatine";
	atp->z = 85;
	atp->M = 210.0;
	atp->radius = 1.4500;
	atp->xsf = 85.0;
	atp->nsf = 0.0;
    }
    // radon
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Rn";
	atp->name = "Radon";
	atp->z = 86;
	atp->M = 222.0;
	atp->radius = 1.4300;
	atp->xsf = 86.0;
	atp->nsf = 0.0;
    }
    // francium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Fr";
	atp->name = "Francium";
	atp->z = 87;
	atp->M = 223.0;
	atp->radius = 2.5000;
	atp->xsf = 87.0;
	atp->nsf = 0.0;
    }
    // radium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ra";
	atp->name = "Radium";
	atp->z = 88;
	atp->M = 226.0;
	atp->radius = 2.1400;
	atp->xsf = 88.0;
	atp->nsf = 10.0;
    }
    // actinium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ac";
	atp->name = "Actinium";
	atp->z = 89;
	atp->M = 227.0;
	atp->radius = 1.8775;
	atp->xsf = 89.0;
	atp->nsf = 0.0;
    }
    // thorium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Th";
	atp->name = "Thorium";
	atp->z = 90;
	atp->M = 232.03811;
	atp->radius = 1.7975;
	atp->xsf = 90.0;
	atp->nsf = 10.31;
    }
    // protactinium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pa";
	atp->name = "Protactinium";
	atp->z = 91;
	atp->M = 231.035882;
	atp->radius = 1.6086;
	atp->xsf = 91.0;
	atp->nsf = 9.13;
    }
    // uranium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "U";
	atp->name = "Uranium";
	atp->z = 92;
	atp->M = 238.028913;
	atp->radius = 1.5683;
	atp->xsf = 92.0;
	atp->nsf = 8.417;
	// isotopes
	// 233-uranium
	itp = new AtomType(*atp);
	itp->symbol = "233-U";
	itp->name = "233-Uranium";
	itp->isotope = true;
	itp->M = 233.039628;
	itp->nsf = 10.1;
	pt_backup.push_back(itp);
	// 234-uranium
	itp = new AtomType(*atp);
	itp->symbol = "234-U";
	itp->name = "234-Uranium";
	itp->isotope = true;
	itp->M = 234.0409456;
	itp->nsf = 12.4;
	pt_backup.push_back(itp);
	// 235-uranium
	itp = new AtomType(*atp);
	itp->symbol = "235-U";
	itp->name = "235-Uranium";
	itp->isotope = true;
	itp->M = 235.0439231;
	itp->nsf = 10.50;
	pt_backup.push_back(itp);
	// 238-uranium
	itp = new AtomType(*atp);
	itp->symbol = "238-U";
	itp->name = "238-Uranium";
	itp->isotope = true;
	itp->M = 238.0507826;
	itp->nsf = 8.407;
	pt_backup.push_back(itp);
    }
    // neptunium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Np";
	atp->name = "Neptunium";
	atp->z = 93;
	atp->M = 237.0;
	atp->radius = 1.0000;
	atp->xsf = 93.0;
	atp->nsf = 10.55;
    }
    // plutonium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Pu";
	atp->name = "Plutonium";
	atp->z = 94;
	atp->M = 244.0;
	atp->radius = 1.0000;
	atp->xsf = 94.0;
	atp->nsf = 7.71;
	// isotopes
	// 240-plutonium
	itp = new AtomType(*atp);
	itp->symbol = "240-Pu";
	itp->name = "240-Plutonium";
	itp->isotope = true;
	itp->M = 240.0538075;
	itp->nsf = 3.5;
	pt_backup.push_back(itp);
	// 242-plutonium
	itp = new AtomType(*atp);
	itp->symbol = "242-Pu";
	itp->name = "242-Plutonium";
	itp->isotope = true;
	itp->M = 242.0587368;
	itp->nsf = 8.1;
	pt_backup.push_back(itp);
    }
    // americium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Am";
	atp->name = "Americium";
	atp->z = 95;
	atp->M = 243.0;
	atp->radius = 1.0000;
	atp->xsf = 95.0;
	atp->nsf = 8.32;
    }
    // curium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cm";
	atp->name = "Curium";
	atp->z = 96;
	atp->M = 247.0;
	atp->radius = 1.0000;
	atp->xsf = 96.0;
	atp->nsf = 9.53;
	// isotopes
	// 246-curium
	itp = new AtomType(*atp);
	itp->symbol = "246-Cm";
	itp->name = "246-Curium";
	itp->isotope = true;
	itp->M = 246.0672176;
	itp->nsf = 9.3;
	pt_backup.push_back(itp);
	// 248-curium
	itp = new AtomType(*atp);
	itp->symbol = "248-Cm";
	itp->name = "248-Curium";
	itp->isotope = true;
	itp->M = 248.072342;
	itp->nsf = 7.7;
	pt_backup.push_back(itp);
    }
    // berkelium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Bk";
	atp->name = "Berkelium";
	atp->z = 97;
	atp->M = 247.0;
	atp->radius = 1.0000;
	atp->xsf = 97.0;
	atp->nsf = 0.0;
    }
    // californium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cf";
	atp->name = "Californium";
	atp->z = 98;
	atp->M = 251.0;
	atp->radius = 1.0000;
	atp->xsf = 98.0;
	atp->nsf = 0.0;
    }
    // einsteinium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Es";
	atp->name = "Einsteinium";
	atp->z = 99;
	atp->M = 252.0;
	atp->radius = 1.0000;
	atp->xsf = 99.0;
	atp->nsf = 0.0;
    }
    // fermium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Fm";
	atp->name = "Fermium";
	atp->z = 100;
	atp->M = 257.0;
	atp->radius = 1.0000;
	atp->xsf = 100.0;
	atp->nsf = 0.0;
    }
    // mendelevium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Md";
	atp->name = "Mendelevium";
	atp->z = 101;
	atp->M = 258.0;
	atp->radius = 1.0000;
	atp->xsf = 101.0;
	atp->nsf = 0.0;
    }
    // nobelium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "No";
	atp->name = "Nobelium";
	atp->z = 102;
	atp->M = 259.0;
	atp->radius = 1.0000;
	atp->xsf = 102.0;
	atp->nsf = 0.0;
    }
    // lawrencium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Lr";
	atp->name = "Lawrencium";
	atp->z = 103;
	atp->M = 262.0;
	atp->radius = 1.0000;
	atp->xsf = 103.0;
	atp->nsf = 0.0;
    }
    // rutherfordium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Rf";
	atp->name = "Rutherfordium";
	atp->z = 104;
	atp->M = 261.0;
	atp->radius = 1.0000;
	atp->xsf = 104.0;
	atp->nsf = 0.0;
    }
    // dubnium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Db";
	atp->name = "Dubnium";
	atp->z = 105;
	atp->M = 262.0;
	atp->radius = 1.0000;
	atp->xsf = 105.0;
	atp->nsf = 0.0;
    }
    // seaborgium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Sg";
	atp->name = "Seaborgium";
	atp->z = 106;
	atp->M = 266.0;
	atp->radius = 1.0000;
	atp->xsf = 106.0;
	atp->nsf = 0.0;
    }
    // bohrium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Bh";
	atp->name = "Bohrium";
	atp->z = 107;
	atp->M = 264.0;
	atp->radius = 1.0000;
	atp->xsf = 107.0;
	atp->nsf = 0.0;
    }
    // hassium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Hs";
	atp->name = "Hassium";
	atp->z = 108;
	atp->M = 277.0;
	atp->radius = 1.0000;
	atp->xsf = 108.0;
	atp->nsf = 0.0;
    }
    // meitnerium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Mt";
	atp->name = "Meitnerium";
	atp->z = 109;
	atp->M = 268.0;
	atp->radius = 1.0000;
	atp->xsf = 109.0;
	atp->nsf = 0.0;
    }
    // darmstadtium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Ds";
	atp->name = "Darmstadtium";
	atp->z = 110;
	atp->M = 281.0;
	atp->radius = 1.0000;
	atp->xsf = 110.0;
	atp->nsf = 0.0;
    }
    // roentgenium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Rg";
	atp->name = "Roentgenium";
	atp->z = 111;
	atp->M = 272.0;
	atp->radius = 1.0000;
	atp->xsf = 111.0;
	atp->nsf = 0.0;
    }
}

// End of file
