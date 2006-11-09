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
* $Id$
*
***********************************************************************/

// Up to date with 1.3.10 Fortran version

#include <fstream>
#include <sstream>
#include <iomanip>

#include "PointsInSphere.h"
#include "PeriodicTable.h"
#include "Atom.h"
#include "StringUtils.h"

#include "pdffit.h"

// Read a number and an eventual comma delimiter or EOF
template<class Type> Type vget(istringstream &fin, char delim)
{
    char c;
    Type val;

    fin >> val;
    //cout << val << " " << c << " " << fin << " " << fin.eof() << endl;

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

template<class Type> Type vget(istringstream &fin)
{
    Type val;

    fin >> val;
    //cout << val << endl;
    return val;
}

double dget(istringstream &fin)
{
    return vget<double>(fin, ',');
}

int iget(istringstream &fin)
{
    return vget<int>(fin, ',');
}

/***********************************************************************
        Read a structure file.
*************************************************************************/
int PdfFit::read_struct(string structfile)
{
    Phase *phase = new Phase;

    try {
        phase->read_struct(nphase+1, structfile);
    }
    catch(Exception e) {
        delete phase;
        // Moved error catching to python bindings.
        throw;
        return 0;
    }

    //phase.off(cr_ipha) = cr_off(cr_ipha-1) + cr_natoms(cr_ipha)*n_at;

    this->phase.push_back(phase);
    total += phase->natoms;
    nphase++;
    setphase(nphase);
    phase->show_lattice();


//  update_cr_dim();
    return 1;
}

/***********************************************************************
        Wed Oct 12 2005 - CLF
        Read a structure from a storage string.
*************************************************************************/
int PdfFit::read_struct_string(char * buffer)
{
    Phase *phase = new Phase;

    try {
        phase->read_struct_string(nphase+1, buffer);
    }
    catch(Exception e) {
        delete phase;
//      cout << "Error in structure"
//          << e.GetMsg() << " --> no new phase allocated\n";
        throw;
        return 0;
    }

    this->phase.push_back(phase);
    total += phase->natoms;
    nphase++;
    setphase(nphase);
    phase->show_lattice();

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
	cout << " Structure file format  : DISCUS (converting B -> Uij) \n";
	Atom::streamformat = Atom::DISCUS;
    }
    else
    {
	cout << " Structure file format  : PDFFIT\n";
	Atom::streamformat = Atom::PDFFIT;
    }

    read_atoms(fstruct);
    // update atom_types
    atom_types.clear();
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
	if (  find(atom_types.begin(), atom_types.end(), ai->atom_type) ==
		atom_types.end() )
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



/************************************************************************************
  compare if string str1 equals string str2, over a length of max(len(str1), minlen)
  Returns 1 if matching, 0 otherwise
*************************************************************************************/
int strcmp(string str1, string str2, int minlen)
{
    int len=str1.size();
    if (len < minlen) return 0;  // string too short to compare
    if (!str2.compare(0,len,str1)) return 1;   // compare if str1 corresponds to substring of str2 over length len
    else return 0;
}

// trims leading spaces
string trim(const string &line)
{
    string trimmed;

    int indx = line.find_first_not_of(" ");
    trimmed = line.substr(indx,line.length()-indx);
    return trimmed;
}

/******************************************************************
    This subroutine reads the header of a structure file
    Wed Oct 12 2005 - CLF
    Changed ifstream to istream to accomodate stringstreams
    as well.
********************************************************************/

void Phase::read_header(istream &fstruct, bool &ldiscus)
{
    //  include     'config.inc'
    //  include     'crystal.inc'
    //  include     'pdf.inc'
    //  include     'errlist.inc'

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

            //------- The maximum number of significant characters depends on the
            //------  length of the character constant befehl.


            //------  Command parameters start at the first character following the blank


            if (befehl == "atoms")
                break;

            //------ Commentary

            else if(befehl[0] == '#')
                continue;

            // ------ Format

            else if(strcmp(befehl,"format",2) )
            {
                string format;

                sline >> format;
                ldiscus = !strcmp(format,"pdffit",3);
            }

            //------ Scale factor (PDFFIT)

            else if(strcmp(befehl,"scale",2) )
            {
                action = "Reading scale factor";
                skal = dget(sline);
                dskal = 0.0;
                if (!sline)
                {
                    throw structureError(action);
                }
            }

            //------    - Peak sharpening (PDFFIT)

            else if(strcmp(befehl,"sharpen",2) )
            {
                action = "reading sharpening parameters";

                // at least 3-parameters must be read without error
                delta2 = dget(sline);
                delta1 = dget(sline);
                srat = dget(sline);

                try { rcut = dget(sline); }

                // if no 4-th parameter available: reshuffle
                catch(vgetException) {
                    rcut = srat;
                    srat = delta1;
                    delta1 = 0;
                }

                ddelta2 = 0.0;
                dsrat = 0.0;
                ddelta1 = 0.0;
            }

            //------ - Space group symbol (only to save it later for DISCUS use)

            else if(strcmp(befehl,"spcgr",1) )
            {
                sline >> spcgr;
            }

            //------    - Title / name for structure

            else if(strcmp(befehl,"title",1) )
            {
                getline(sline, name);  // does not trim leading spaces
                name = trim(name);
            }

            //------    - Cell constants

            else if(strcmp(befehl,"cell",1) )
            {
                action = "reading unit cell parameters";
                a0[0] = dget(sline);
                a0[1] = dget(sline);
                a0[2] = dget(sline);
                win[0] = dget(sline);
                win[1] = dget(sline);
                win[2] = dget(sline);
            }

            //------    - Standard deviation of cell constants

            else if(strcmp(befehl,"dcell",2) )
            {
                action = "reading standard deviation of unit cell parameters";
                da0[0] = dget(sline);
                da0[1] = dget(sline);
                da0[2] = dget(sline);
                dwin[0] = dget(sline);
                dwin[1] = dget(sline);
                dwin[2] = dget(sline);
            }

            //------    - Crystal dimensions and number of atoms per unit cell 'ncell'

            else if(strcmp(befehl,"ncell",1) )
            {
                action = "reading # atoms/unit cell";
                icc[0] = iget(sline);
                icc[1] = iget(sline);
                icc[2] = iget(sline);
                ncatoms = iget(sline);
            }

            else
            {
                cout << " ****WARN**** Unknown keyword: " << befehl << " (ignored) ****\n";
            }
        } // end try
        // catch vget-exception and throw the specific exception
        catch(vgetException e) { throw structureError(action+e.GetMsg()); }
    }
}

void Phase::read_atoms(istream& fstruct)
{
    Atom a;
    while (fstruct >> a)
    {
        this->atom.push_back(a);
        natoms++;
    }
    return;
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
        //warning("save_struct: phase does not exist\n");
        throw unassignedError("phase does not exist");
    }
    else
    {
        bool ldiscus=false;
        ofstream fout;

        if (!strucfile.empty())
        {
            if (ldiscus)
                cout << " Saving structure (DISCUS format) phase " << ip
            << " to file : " << strucfile << endl;
            else
                cout << " Saving structure phase " << ip << " to file : "
            << strucfile << endl;

            phase[ip-1]->save_struct(outfilestream);

            fout.open(strucfile.c_str());
            if (!fout) {
                //warning("save_struct: cannot create output file");
                throw IOError("cannot create output file");
                return "";
            }
            fout << outfilestream.str();
            fout.close();
        }
        else
        {
            phase[ip-1]->save_struct(outfilestream);
//            cout << outfilestream.str();
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

    if(!ldis)
    {
	fout << "format pdffit" << endl;
	fout << "scale  " << setw(9) << skal << endl;
	fout << "sharp  " << setw(9) << delta2 << ", " << setw(9) << delta1 << ", "
	    << setw(9) << srat << ", " << setw(9) << rcut << endl;
    }

    fout << "spcgr  " << spcgr << endl;

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
	    fout << setw(4) << left << toupper(ai->atom_type->symbol);
	    fout << right << setprecision(8);
	    for (int i=0; i<3; i++) fout << setw(18) << ai->pos[i];
	    fout << setw(13) << dw << endl;
	}
	else
	{
	    fout << setw(4) << left << toupper(ai->atom_type->symbol);
	    fout << right << setprecision(8);
	    for (int i=0; i<3; i++) fout << setw(18) << ai->pos[i];
	    fout << setw(13) << setprecision(4) << ai->occ << endl;

	    fout << "    ";
	    fout << setprecision(8);
	    for (int i=0; i<3; i++) fout << setw(18) << ai->dpos[i];
	    fout << setw(13) << setprecision(4) << ai->docc << endl;

	    fout << "    ";
	    fout << setprecision(8);
	    for (int i=0; i<3; i++) fout << setw(18) << ai->u[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=0; i<3; i++) fout << setw(18) << ai->du[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=3; i<6; i++) fout << setw(18) << ai->u[i];
	    fout << endl;

	    fout << "    ";
	    for (int i=3; i<6; i++) fout << setw(18) << ai->du[i];
	    fout << endl;
	}
    }
}


/***************************************
c   Calculate bond angles with errors
****************************************/
double PdfFit::bond_angle(int ia, int ja, int ka)
{
    if( !curphase )
    {
        throw unassignedError("Must read structure first");
        return 0;
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
	xyz[i] = remainder(xyz[i], icc[i]);
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

double Phase::bond_angle(int ia, int ja, int ka)
{
    double x[3], y[3], dx[3], dy[3], xx, yy, xy, dxx, dyy, dxy, arg, darg, ang, dang;


    if ( (ia < 1) || (ia > natoms) || (ja < 1) || (ja > natoms)
        || (ka < 1) || (ka > natoms))
    {
        //cout << "Incorrect atom number\n";
        stringstream eout;
        eout << "Incorrect atom number(s): " << ia << ", " << ja << ", " << ka;
        throw ValueError(eout.str());
        return 0;
    }
    if ( ia == ja || ia == ka || ja == ka )
    {
        stringstream eout;
        eout << "All atoms must be different: " << ia << ", ";
        eout << ja << ", " << ka;
        throw ValueError(eout.str());
        return 0;
    }

    ia--; ja--; ka--;

    for (int i=0; i<3; i++)
    {
        x[i] = atom[ja].pos[i] - atom[ia].pos[i];
        y[i] = atom[ja].pos[i] - atom[ka].pos[i];
        dx[i] = atom[ja].dpos[i] + atom[ia].dpos[i];
        dy[i] = atom[ja].dpos[i] + atom[ka].dpos[i];
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
        dang = abs(1.0/sqrt(1.0-arg*arg)/rad*darg);
    else
        dang = 0.0;

    cout << "   "
       << toupper(atom[ia].atom_type->symbol) << " (#" << ia+1 << ")"
       <<" - "
       << toupper(atom[ja].atom_type->symbol) << " (#" << ja+1 << ")"
       <<" - "
       << toupper(atom[ka].atom_type->symbol) << " (#" << ka+1 << ")"
       << "   =   "  << putxdx(ang,dang) << " degrees\n";

   return ang;
}

/***************************************
c   Calculate bond lengths with errors
****************************************/
double PdfFit::bond_length_atoms(int ia, int ja)
{
    if( !curphase )
    {
        throw unassignedError("Must read structure first");
        return 0;
    }
    return curphase->bond_length_atoms(ia, ja);
}

double Phase::bond_length_atoms(int ia, int ja)
{
    double d[3], dd[3], dist, ddist;

    // -- Simple distance between given atoms

    if ( (ia < 1) || (ia > natoms) || (ja < 1) || (ja > natoms) )
    {

	stringstream eout;
	eout << "Incorrect atom number(s): " << ia << ", " << ja;
	throw ValueError(eout.str());
	throw ValueError("Incorrect atom number");
	return 0;
    }

    ia--; ja--;

    for (int jj=0; jj<3; jj++)
    {
	d[jj] = atom[ia].pos[jj] - atom[ja].pos[jj];
	dd[jj] = atom[ia].dpos[jj] + atom[ja].dpos[jj];
    }
    make_nearest(d);
    dist = sqrt(skalpro(d,d));
    ddist = 0.5/dist * dskalpro(d,d,dd,dd);
    cout << "   " << toupper(atom[ia].atom_type->symbol) << " (#" << ia+1 << ")"
	<< " - " << toupper(atom[ja].atom_type->symbol) << " (#" << ja+1 <<
	")   =   " << putxdx(dist,ddist) << " A" << endl;
    return dist;
}


void PdfFit::bond_length_types(const string& symi, const string& symj,
	double bmin, double bmax)
{
    if(!curphase)
    {
        throw unassignedError("Must read structure first");
    }
    curphase->bond_length_types(symi, symj, bmin, bmax);
}

void Phase::bond_length_types(string symi, string symj,
	double bmin, double bmax)
{
    bool lfound;
    double d[3], dd[3], dist, ddist;
    set<size_t> iselection, jselection;
    iselection = selectAtomsOf(symi);
    jselection = selectAtomsOf(symj);

    lfound = false;

    // ---- Get all bonds in specified range

    cout << "(" << toupper(symi);
    cout << "," << toupper(symj) << ")";
    cout << " bond lengths in [" << bmin << "A," << bmax << "A]";
    cout << " for current phase : \n";

    // calculate range for PointsInSphere sequencer
    // (negative rsphmin is no problem)
    double buffzone = circum_diameter();
    double rsphmin = bmin - buffzone;
    double rsphmax = bmax + buffzone;
    PointsInSphere sph( rsphmin, rsphmax, a0[0]*icc[0],
	    a0[1]*icc[1], a0[2]*icc[2],
	    win[0], win[1], win[2] );

    // -- Loop over all atoms within the crystal

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
		    ddist = 0.5/dist * dskalpro(d,d,dd,dd);
		    cout << "   " << toupper(atom[*ia].atom_type->symbol) <<
			" (#" << *ia+1 << ")" << " - " <<
			toupper(atom[*ja].atom_type->symbol) <<
			" (#" << *ja+1 << ")   =   " <<
			putxdx(dist,ddist) << " A" << endl;
		    lfound = true;
		}
	    }
	}
    }

    cout << endl;

    if (!lfound) cout << "   *** No pairs found ***\n";
}

set<size_t> Phase::selectAtomsOf(string symbol)
{
    set<size_t> selection;
    symbol = toupper(symbol);
    if (symbol == "ALL")
    {
	for (size_t i = 0; i != size_t(natoms); ++i)  selection.insert(i);
	return selection;
    }
    // here we need to find AtomType
    PeriodicTable* pt = PeriodicTable::instance();
    AtomType* atp;
    try
    {
	atp = pt->lookup(symbol);
    }
    catch (runtime_error e)
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
