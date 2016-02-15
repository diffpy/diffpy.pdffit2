/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* PdfFit and Phase methods for reading and saving the structure,
* and for calculating bond lengths and angles.
*
* Comments:
*
***********************************************************************/

// Up to date with 1.3.10 Fortran version

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "PointsInSphere.h"
#include "LocalPeriodicTable.h"
#include "Atom.h"
#include "StringUtils.h"
#include "PairDistance.h"
#include "MathUtils.h"
#include "pdffit.h"

using NS_PDFFIT2::pout;


/***********************************************************************
* local helper routines
***********************************************************************/

namespace {

// Read a number and an eventual comma delimiter or EOF
template<class Type> Type vget(istringstream &fin, char delim)
{
    char c;
    Type val;

    fin >> val;

    // Return if reading error
    if (!fin) {
        string line;
        fin.clear();
        fin >> line;
        throw vgetException(line);
    }
    if (!fin) return 0;

    // check on end of file before looking for delimiter
    if (fin.eof()) return val;

    // read eventual delimiter
    fin >> c;

    // if a character c has been read and it is not the expected
    // delimiter: put back for next reading
    if (fin && c != delim)  fin.unget();

    // clear any error in case no delimiter could be read
    if (!fin)	fin.clear();

    return val;
}

// read space delimited value
template<class Type> Type vget(istringstream &fin)
{
    Type val;
    fin >> val;
    return val;
}

// read space or comma delimited double
double dget(istringstream &fin)
{
    return vget<double>(fin, ',');
}

// read space or comma delimited integer
int iget(istringstream &fin)
{
    return vget<int>(fin, ',');
}


// strip leading spaces
string lstrip(const string &line)
{
    string naked;
    string::size_type i = line.find_first_not_of(" \t");
    if (i != string::npos)  naked = line.substr(i);
    return naked;
}

// substitute all occurences of literal pattern with new string
void substitute(string& s, const string& pat, const string& sub)
{
    string::size_type p;
    string::size_type start = 0;
    for (p = s.find(pat, start); p != string::npos; p = s.find(pat, start))
    {
        s = s.replace(p, pat.size(), sub);
        start = p + sub.size();
    }
}


}   // local namespace


/***********************************************************************
* Read a structure file.
***********************************************************************/
int PdfFit::read_struct(string structfile)
{
    Phase* ph = new Phase();
    try {
        ph->read_struct(nphase+1, structfile);
    }
    catch(Exception e) {
        delete ph;
        // Moved error catching to python bindings.
        throw;
    }
    this->phase.push_back(ph);
    this->total += ph->natoms;
    this->nphase++;
    this->selphaseForEachDataSet(ph);
    this->setphase(this->nphase);
    ph->show_lattice();
    return 1;
}

/***********************************************************************
        Wed Oct 12 2005 - CLF
        Read a structure from a storage string.
*************************************************************************/
int PdfFit::read_struct_string(char * buffer)
{
    Phase* ph = new Phase();
    try {
        ph->read_struct_string(nphase+1, buffer);
    }
    catch(Exception e) {
        delete ph;
        throw;
    }
    this->phase.push_back(ph);
    this->total += ph->natoms;
    this->nphase++;
    this->selphaseForEachDataSet(ph);
    this->setphase(this->nphase);
    ph->show_lattice();
    return 1;
}

void Phase::read_struct(int _iphase, string structfile)
{
    ifstream fstruct;

    fstruct.open(structfile.c_str());
    if (!fstruct) throw IOError("File does not exist");

    read_struct_stream(_iphase, fstruct);
}

void Phase::read_struct_string(int _iphase, char * buffer)
{
    istringstream fstruct(buffer);
    read_struct_stream(_iphase, fstruct);
}

void Phase::read_struct_stream(int _iphase, istream& fstruct)
{
    double tot;
    bool ldiscus;

    iphase = _iphase;
    natoms = 0;

    read_header(fstruct, ldiscus);

    if (ldiscus)
    {
	*pout << " Structure file format  : DISCUS (converting B -> Uij) \n";
	Atom::streamformat = Atom::DISCUS;
    }
    else
    {
	*pout << " Structure file format  : PDFFIT\n";
	Atom::streamformat = Atom::PDFFIT;
    }

    read_atoms(fstruct);
    // update atom_types
    atom_types.clear();
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	if (!count(atom_types.begin(), atom_types.end(), ai->atom_type))
	{
	    atom_types.push_back(ai->atom_type);
	}
    }

    lattice();

    tot = icc[0]*icc[1]*icc[2]*ncatoms;

    if (tot != natoms)
    {
	throw structureError("Inconsistent # of atoms in structure");
    }
}


/******************************************************************
    This subroutine reads the header of a structure file
    Wed Oct 12 2005 - CLF
    Changed ifstream to istream to accomodate stringstreams
    as well.
********************************************************************/

void Phase::read_header(istream &fstruct, bool &ldiscus)
{
    string ier_msg;
    string befehl, line;


    // initialize format type to discuss format
    ldiscus = true;

    /*  parse structure file header and fill phase variables */

    while (getline(fstruct, line))
    {
        istringstream sline(line);
        string action;

        // try to read a command from structure file
        // if an error occurs, an exception will be caught
        try{
            sline >> befehl;

	    // get out if we get to atom positions
            if (befehl == "atoms")	break;

            // skip comments, i.e., when befehl starts with '#'
            else if (befehl.find('#') == 0)  continue;

            // format
            else if (befehl == "format")
            {
                string format;
                sline >> format;
                ldiscus = (format != "pdffit");
            }

            // scale factor (PDFFIT)
            else if (befehl == "scale")
            {
                action = "Reading scale factor";
                this->pscale = dget(sline);
                this->dpscale = 0.0;
                if (!sline)
                {
                    throw structureError(action);
                }
            }

            // peak sharpening factors (PDFFIT)
            else if (befehl == "sharp")
            {
                action = "reading sharpening parameters";
		double v0 = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
                // at least 3-parameters must be read without error
		v0 = dget(sline);
		v1 = dget(sline);
		v2 = dget(sline);
		// we have new format if we can read the 4th parameter
		try {
		    v3 = dget(sline);
		    delta2 = v0;
		    delta1 = v1;
		    sratio = v2;
		    rcut = v3;
		}
		// if reading of 4th parameter fails, assume old format
                catch(vgetException) {
		    delta2 = v0;
		    delta1 = 0.0;
		    sratio = v1;
		    rcut = v2;
		}
                ddelta2 = 0.0;
                dsratio = 0.0;
                ddelta1 = 0.0;
            }

            // space group symbol (only to save it later for DISCUS use)
            else if (befehl == "spcgr")
            {
                sline >> spcgr;
            }

            // particle shape corrections
            else if (befehl == "shape")
            {
                action = "reading particle shape correction data";
                string shapedata;
                getline(sline, shapedata);
                substitute(shapedata, ",", " ");
                istringstream shapestream(shapedata);
                string w;
                shapestream >> w;
                if (w == "sphere")
                {
                    this->spdiameter = dget(shapestream);
                }
                else if (w == "stepcut")
                {
                    this->stepcut = dget(shapestream);
                }
                else
                {
                    ostringstream emsg;
                    emsg << " Unknown shape correction: " << w;
                    throw structureError(emsg.str());
                }
            }

            // title / name for structure
            else if (befehl == "title")
            {
                getline(sline, name);
		// getline keeps leading whitespace
                name = lstrip(name);
            }

            // cell constants
            else if (befehl == "cell")
            {
                action = "reading unit cell parameters";
                a0[0] = dget(sline);
                a0[1] = dget(sline);
                a0[2] = dget(sline);
                win[0] = dget(sline);
                win[1] = dget(sline);
                win[2] = dget(sline);
            }

            // standard deviation of cell constants
            else if (befehl == "dcell")
            {
                action = "reading standard deviation of unit cell parameters";
                da0[0] = dget(sline);
                da0[1] = dget(sline);
                da0[2] = dget(sline);
                dwin[0] = dget(sline);
                dwin[1] = dget(sline);
                dwin[2] = dget(sline);
            }

            // crystal dimensions and number of atoms per unit cell 'ncell'
            else if (befehl == "ncell")
            {
                action = "reading # atoms/unit cell";
                icc[0] = iget(sline);
                icc[1] = iget(sline);
                icc[2] = iget(sline);
                ncatoms = iget(sline);
            }

	    // show warning message otherwise
            else
            {
                *pout << " ****WARN**** Unknown keyword: " <<
                    befehl << " (ignored) ****\n";
            }
        }   // end of try
        // catch vget-exception and throw the specific exception
        catch(vgetException e) { throw structureError(action+e.GetMsg()); }
    }
}

void Phase::read_atoms(istream& fstruct)
{
    Atom a;
    while (fstruct >> a)
    {
        reassign_atom_type(&a);
        this->atom.push_back(a);
        natoms++;
    }
    return;
}


void Phase::reassign_atom_type(Atom* ap)
{
    LocalPeriodicTable* lpt = getPeriodicTable();
    const string& smbl = ap->atom_type->symbol;
    ap->atom_type = lpt->symbol(smbl);
}



/*******************************************
c------ - Save structure for given phase
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the
    saved file. Actually saving the file is
    optional.
*******************************************/
string PdfFit::show_struct(int ip)
{
    string filestring = save_struct(ip,"");
    return filestring;
}

string PdfFit::save_struct(int ip, string strucfile)
{

    stringstream outfilestream;

    if ( (ip < 1) || (ip > nphase) )
    {
        throw unassignedError("phase does not exist");
    }
    else
    {
        bool ldiscus=false;
        ofstream fout;

        if (!strucfile.empty())
        {
            if (ldiscus)
                *pout << " Saving structure (DISCUS format) phase " << ip
            << " to file : " << strucfile << endl;
            else
                *pout << " Saving structure phase " << ip << " to file : "
            << strucfile << endl;

            phase[ip-1]->save_struct(outfilestream);

            fout.open(strucfile.c_str());
            if (!fout) {
                throw IOError("cannot create output file");
            }
            fout << outfilestream.str();
            fout.close();
        }
        else
        {
            phase[ip-1]->save_struct(outfilestream);
        }
    }

    return outfilestream.str();
}

/***************************************************************
c   This subroutine saves the structure and/or the unit cell
c   onto a file. The format uses keyword description.
****************************************************************/
template <class Stream> void Phase::save_struct(Stream &fout)
{
    const double fac = 8.0*M_PI*M_PI /3.0;
    bool ldis=false;

    // -- Write new type of structure file

    fout << "title  " << name << endl;

    fout << setprecision(6) << fixed;

    if (!ldis)
    {
	fout << "format pdffit" << endl;
	fout << "scale  " << setw(9) << pscale << endl;
	fout << "sharp  " << setw(9) << delta2 << ", " << setw(9) << delta1 << ", "
	    << setw(9) << sratio << ", " << setw(9) << rcut << endl;
    }

    fout << "spcgr   " << spcgr << endl;

    if (spdiameter > 0.0)
    {
        fout << "shape   sphere, " << spdiameter << endl;
    }

    if (stepcut > 0.0)
    {
        fout << "shape   stepcut, " << stepcut << endl;
    }

    fout << "cell   ";
    for (int i=0; i<3; i++)
	fout << setw(9) << a0[i] << ", ";

    for (int i=0; i<3; i++)
    {
	fout << setw(9) << win[i];
	if (i!=2)
	    fout << ", ";
	else
	    fout << endl;
    }

    if (!ldis)
    {
	fout << "dcell  ";
	for (int i=0; i<3; i++)
	    fout << setw(9) << da0[i] << ", ";

	for (int i=0; i<3; i++)
	{
	    fout << setw(9) << dwin[i];
	    if (i!=2)
		fout << ", ";
	    else
		fout << endl;
	}
    }

    fout << "ncell  ";
    for (int i=0; i<3; i++)
	fout << setw(9) << icc[i] << ", ";
    fout << setw(9) << ncatoms << endl;

    fout << "atoms" << endl;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {

	if (ldis)
	{
	    double dw = fac*(ai->u[1]+ai->u[2]+ai->u[3]);
	    fout << setw(4) << left << ai->atom_type->symbol;
	    fout << right << setprecision(8);
	    for (int i=0; i<3; i++) fout << ' ' << setw(17) << ai->pos[i];
	    fout << ' ' << setw(12) << dw << endl;
	}
	else
	{
	    fout << setw(4) << left << ai->atom_type->symbol;
	    fout << right << setprecision(8);
	    for (int i=0; i<3; i++) fout << ' ' << setw(17) << ai->pos[i];
	    fout << ' ' << setw(12) << setprecision(4) << ai->occ << endl;

	    fout << "    ";
	    fout << setprecision(8);
	    for (int i=0; i<3; i++) fout << ' ' << setw(17) << ai->dpos[i];
	    fout << ' ' << setw(12) << setprecision(4) << ai->docc << endl;

	    fout << "    ";
	    fout << setprecision(8);
	    for (int i=0; i<3; i++) fout << ' ' << setw(17) << ai->u[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=0; i<3; i++) fout << ' ' << setw(17) << ai->du[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=3; i<6; i++) fout << ' ' << setw(17) << ai->u[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=3; i<6; i++) fout << ' ' << setw(17) << ai->du[i];
	    fout << endl;
	}
    }
}


/***************************************
c   Calculate bond angles with errors
****************************************/
pair<double,double> PdfFit::bond_angle(int ia, int ja, int ka)
{
    if (!curphase)
    {
        throw unassignedError("Must read structure first");
    }
    return curphase->bond_angle(ia, ja, ka);
}

/***********************************************************************
* shift to equivalent lattice position that is nearest to the origin
***********************************************************************/
void Phase::make_nearest(double xyz[3])
{
    // first shift to the first unit cell
    for (int i = 0; i !=3; ++i)
    {
	// MS compatibility - use fmod instead of remainder
	xyz[i] = fmod(xyz[i], icc[i]);
    }
    // that is all in orthogonal cell or if we get to the origin
    if ( (xyz[0] == 0.0 && xyz[1] == 0.0 && xyz[2] == 0.0) ||
	 (win[0] == 90.0 && win[1] == 90.0 && win[2] == 90.0) )
    {
	return;
    }
    // otherwise we need to check all cells around
    // first, let us shift to octant with xyz[i] <= 0.0
    for (int i = 0; i != 3; ++i)
    {
	if (xyz[i] > 0.0)   xyz[i] -= icc[i];
    }
    double nearest[3] = { xyz[0], xyz[1], xyz[2] };
    double mindsquare = skalpro(nearest, nearest);
    double test[3];
    for (test[0] = xyz[0]; test[0] < icc[0]; test[0] += icc[0])
    {
	for (test[1] = xyz[1]; test[1] < icc[1]; test[1] += icc[1])
	{
	    for (test[2] = xyz[2]; test[2] < icc[2]; test[2] += icc[2])
	    {
		double dsquare = skalpro(test,test);
		if (dsquare < mindsquare)
		{
		    copy(test, test+3, nearest);
		    dsquare = mindsquare;
		}
	    }
	}
    }
    copy(nearest, nearest+3, xyz);
}

pair<double,double> Phase::bond_angle(int ia, int ja, int ka)
{
    double x[3], y[3], dx[3], dy[3], xx, yy, xy, dxx, dyy, dxy, arg, darg, ang, dang;


    if ( (ia < 1) || (ia > natoms) || (ja < 1) || (ja > natoms)
        || (ka < 1) || (ka > natoms))
    {
        stringstream eout;
        eout << "Incorrect atom number(s): " << ia << ", " << ja << ", " << ka;
        throw ValueError(eout.str());
    }
    if ( ia == ja || ia == ka || ja == ka )
    {
        stringstream eout;
        eout << "All atoms must be different: " << ia << ", ";
        eout << ja << ", " << ka;
        throw ValueError(eout.str());
    }

    Atom& ai = atom[ia - 1];
    Atom& aj = atom[ja - 1];
    Atom& ak = atom[ka - 1];

    for (int i=0; i<3; i++)
    {
        x[i] = aj.pos[i] - ai.pos[i];
        y[i] = aj.pos[i] - ak.pos[i];
        dx[i] = aj.dpos[i] + ai.dpos[i];
        dy[i] = aj.dpos[i] + ak.dpos[i];
    }

    make_nearest(x);
    make_nearest(y);
    xx = sqrt(skalpro(x,x));
    yy = sqrt(skalpro(y,y));
    xy = skalpro(x,y);

    dxx = 0.5/xx*dskalpro(x,x,dx,dx);
    dyy = 0.5/yy*dskalpro(y,y,dy,dy);
    dxy = dskalpro(x,y,dx,dy);
    arg = xy/(xx*yy);
    ang = acosd(arg);
    darg = (1.0/(xx*yy)*dxy + arg/xx*dxx +arg/yy*dyy);
    if (arg != 1.0)
        dang = fabs(1.0/sqrt(1.0-arg*arg)/rad*darg);
    else
        dang = 0.0;

    pair<double,double> rv(ang, dang);

    return rv;
}

/***************************************
c   Calculate bond lengths with errors
****************************************/
PairDistance PdfFit::bond_length_atoms(int ia, int ja)
{
    if (!curphase)
    {
        throw unassignedError("Must read structure first");
    }
    return curphase->bond_length_atoms(ia, ja);
}

PairDistance Phase::bond_length_atoms(int ia, int ja)
{
    double d[3], dd[3], dist, ddist;

    // -- Simple distance between given atoms

    if ( (ia < 1) || (ia > natoms) || (ja < 1) || (ja > natoms) )
    {

	stringstream eout;
	eout << "Incorrect atom number(s): " << ia << ", " << ja;
	throw ValueError(eout.str());
    }

    Atom& ai = atom[ia-1];
    Atom& aj = atom[ja-1];

    for (int jj=0; jj<3; jj++)
    {
	d[jj] = ai.pos[jj] - aj.pos[jj];
	dd[jj] = ai.dpos[jj] + aj.dpos[jj];
    }
    make_nearest(d);
    dist = sqrt(skalpro(d,d));
    ddist = (dist > 0) ? 0.5/dist * dskalpro(d,d,dd,dd) : 0.0;

    PairDistance pd;
    pd.dij = dist;
    pd.ddij = ddist;
    pd.i = ia;
    pd.j = ja;

    return pd;
}


vector<PairDistance> PdfFit::bond_length_types(string symi, string symj,
	double bmin, double bmax)
{
    if (!curphase)
    {
        throw unassignedError("Must read structure first");
    }
    return curphase->bond_length_types(symi, symj, bmin, bmax);
}

vector<PairDistance> Phase::bond_length_types(string symi, string symj,
	double bmin, double bmax)
{
    double d[3], dd[3], dist, ddist;
    set<size_t> iselection, jselection;
    iselection = selectAtomsOf(symi);
    jselection = selectAtomsOf(symj);

    // ---- Get all bonds in specified range

    // calculate range for PointsInSphere sequencer
    // (negative rsphmin is no problem)
    double buffzone = circum_diameter();
    double rsphmin = bmin - buffzone;
    double rsphmax = bmax + buffzone;
    PointsInSphere sph( rsphmin, rsphmax, a0[0]*icc[0],
	    a0[1]*icc[1], a0[2]*icc[2],
	    win[0], win[1], win[2] );

    // -- Loop over all atoms within the crystal
    vector<PairDistance> rv;
    set<size_t>::iterator ia, ja;
    for (ia = iselection.begin(); ia != iselection.end(); ++ia)
    {
	for (ja = jselection.begin(); ja != jselection.end(); ++ja)
	{
	    for (sph.rewind(); !sph.finished(); sph.next())
	    {
		for (int jj=0; jj<3; jj++)
		{
		    d[jj] = atom[*ia].pos[jj] - atom[*ja].pos[jj] -
			    sph.mno[jj]*icc[jj];
		    dd[jj] = atom[*ia].dpos[jj] + atom[*ja].dpos[jj];
		}
		dist = sqrt(skalpro(d,d));
		if ( (dist >= bmin) && (dist <= bmax) )
		{
		    ddist = (dist > 0) ? 0.5/dist * dskalpro(d,d,dd,dd) : 0.0;
		    PairDistance pd;
		    pd.dij = dist;
		    pd.ddij = ddist;
		    pd.i = *ia + 1;
		    pd.j = *ja + 1;
		    rv.push_back(pd);
		}
	    }
	}
    }
    stable_sort(rv.begin(), rv.end());
    return rv;
}

set<size_t> Phase::selectAtomsOf(string symbol)
{
    set<size_t> selection;
    if (toupper(symbol) == "ALL")
    {
	for (size_t i = 0; i != size_t(natoms); ++i)  selection.insert(i);
	return selection;
    }
    // here we need to find AtomType
    LocalPeriodicTable* lpt = getPeriodicTable();
    const AtomType* atp;
    try
    {
	atp = lpt->lookup(symbol);
    }
    catch (ValueError e)
    {
	ostringstream emsg;
	emsg << "Incorrect atom type '" << symbol << "'";
	throw ValueError(emsg.str());
    }
    for (size_t i = 0; i != size_t(natoms); ++i)
    {
	if (atom[i].atom_type == atp)	selection.insert(i);
    }
    return selection;
}

// End of file
