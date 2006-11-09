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
* $Id$
*
***********************************************************************/

#include <sstream>
#include <stdexcept>

#include "PeriodicTable.h"

using namespace std;

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

AtomType* PeriodicTable::atomicNumber(size_t z)
{
    if (z == 0 || z > max_z)
    {
	ostringstream emsg;
	emsg << "Element number " << z << " is not defined.";
	throw runtime_error(emsg.str());
    }
    return pt_public[z - 1];
}

AtomType* PeriodicTable::name(const string& s)
{
    map<string,AtomType*>::iterator ii;
    ii = name_index.find(s);
    if (ii == name_index.end())
    {
	ostringstream emsg;
	emsg << "Element '" << s << "' is not defined.";
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
	emsg << "Element '" << s << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    return ii->second;
}

AtomType* PeriodicTable::lookup(string s)
{
    // force standard case
    if (!s.empty())
    {
	string::iterator sii = s.begin();
	*sii = toupper(*sii);
	for (sii++; sii != s.end(); ++sii)  *sii = tolower(*sii);
    }
    map<string,AtomType*>::iterator ii;
    ii = symbol_index.find(s);
    if (    ii == symbol_index.end() &&
	    (ii = name_index.find(s)) == name_index.end() )
    {
	ostringstream emsg;
	emsg << "Element '" << s << "' is not defined.";
	throw runtime_error(emsg.str());
    }
    return ii->second;
}

void PeriodicTable::defAtomType(const AtomType atp)
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
    // Refs:
    // 1. Albert-Jose Dianoux, Gerry Lander, Neutron Data Booklet,
    //    Second Edition, ILL 2003
    // 2. ionic radii: http://www.fhi-berlin.mpg.de/th/balsac/balm.47.html
    AtomType* atp;
    // hydrogen
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "H";
	atp->name = "Hydrogen";
	atp->z = 1;
	atp->M = 1.007947;
	atp->radius = 0.4350;
	atp->xsf = 1.0;
	atp->nsf = -3.740911;
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
    // caesium
    atp = *pt_backup.insert(pt_backup.end(), new AtomType());
    {
	atp->symbol = "Cs";
	atp->name = "Caesium";
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
    // number of predefined elements
    max_z = pt_backup.size();
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

// End of file
