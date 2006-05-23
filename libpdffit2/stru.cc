// Up to date with 1.3.10 Fortran version

#include <fstream>
#include <sstream>
#include <iomanip>
#include "math.h"
#include "pdffit.h"
#include "stdarg.h"
#include "PointsInSphere.h"


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
    if ( fin && (c != delim) ) fin.unget();

    // //clear any error in case no delimiter could be read before EOF
    // fin.clear no longer necessary since we introduced exception in vget
    //if (fin.eof()) fin.clear();

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



/**************************
    Main read command
**************************/
/*void PdfFit::do_read(string type, double lpara[4])
{

    // Possible read-types: data, diff


// -----    Reading a PDF difference file

      else if(type == "diff")
      {
          read_data(true);
          fit_update();
      }
}*/



/***********************************************************************
        Read a structure file.
*************************************************************************/
// WHY IS LATTICE CALLED TWICE??
int PdfFit::read_struct(string structfile)
{
    Phase *phase = new Phase;

    //call do_build_name (ianz,cpara,lpara,werte,maxw,1)  possibility to number
    //the file

    try {
        phase->read_struct(nphase+1, structfile, total);
    }
    catch(Exception e) {
        delete phase;
        // Moved error catching to python bindings.
//      cout << "Error in structure file <" << structfile << ">: "
//          << e.GetMsg() << " --> no new phase allocated\n";
        throw;
        return 0;
    }

    //phase.off(cr_ipha) = cr_off(cr_ipha-1) + cr_natoms(cr_ipha)*n_at;

    this->phase.push_back(phase);
    nphase++;
    setphase(nphase);
    phase->lattice(true);


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
        phase->read_struct_string(nphase+1, buffer, total);
    }
    catch(Exception e) {
        delete phase;
//      cout << "Error in structure"
//          << e.GetMsg() << " --> no new phase allocated\n";
        throw;
        return 0;
    }

    this->phase.push_back(phase);
    nphase++;
    setphase(nphase);
    phase->lattice(true);

    return 1;
}

/*
    Wed Oct 12 2005 - CLF
    Moved bulk of functionality to read_struct_string.
*/

void Phase::read_struct(int _iphase, string structfile, int &total)
{
    ifstream fstruct;

    fstruct.open(structfile.c_str());
    if (!fstruct) throw IOError("File does not exist");

    // Read the file into a buffer and send it off to read_struct_string function
    char * buffer;
    int length;
    fstruct.seekg( 0, ios::end );
    length = fstruct.tellg();
    fstruct.seekg( 0, ios::beg );
    buffer = new char [length];
    fstruct.read( buffer, length );
    fstruct.close();
    read_struct_string( _iphase, buffer, total );

}

/*
    Wed Oct 12 2005 - CLF
    Added the ability to read a structure stored in a c-style
    string.
*/

void Phase::read_struct_string(int _iphase, char * buffer, int &total)
{
    istringstream fstruct( buffer );
    double tot;
    bool ldiscus;

    iphase = _iphase;
    natoms = 0;

    read_header(fstruct, ldiscus);

    lattice(false);

    if (ldiscus)
      cout << " Structure file format  : DISCUS (converting B -> Uij) \n";
    else
      cout << " Structure file format  : PDFFIT\n";

    read_atoms(fstruct,ldiscus);

    tot = icc[0]*icc[1]*icc[2]*ncatoms;

    if (tot != natoms) throw structureError("Inconsistent # of atoms in structure");

    total += natoms;

    //fstruct.close();
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
                delta = dget(sline);
                gamma = dget(sline);
                srat = dget(sline);

                try { rcut = dget(sline); }

                // if no 4-th parameter available: reshuffle
                catch(vgetException) {
                    rcut = srat;
                    srat = gamma;
                    gamma = 0;
                }

                ddelta = 0.0;
                dsrat = 0.0;
                dgamma = 0.0;
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

/*********************************************************
    Wed Oct 12 2005 - CLF
    Changed ifstream to istream to accomodate stringstreams
    as well.
**********************************************************/
void read3(istream &fin, double *val)
{
    fin >> val[0] >> val[1] >> val[2];
}

/*
    Wed Oct 12 2005 - CLF
    Changed ifstream to istream to accomodate stringstreams
    as well.
*/
void Phase::read_atoms(istream &fstruct, bool ldiscus)
{

    //include   'config.inc'
    //include   'crystal.inc'
    //include   'wink.inc'
    //include   'errlist.inc'

    Atom atom;

    while (!atom.read_atom(fstruct, ldiscus, this))
    {
        this->atom.push_back(atom);
        natoms++;
    }
    // Wed Oct 12 2005 - CLF
    // This function should not close any streams
    //fstruct.close();
    return;
}

/***********************************************************
    Reads the list of atoms into the crystal array
    Wed Oct 12 2005 - CLF
    Changed ifstream to istream to accomodate stringstreams
    as well.
************************************************************/
int Atom::read_atom(istream &fstruct, bool ldiscus, void *_phase)
{
    Phase *phase;
    double x, y, z, dw;
    string type;
    int j;
    const double fac = 1.0 / (8.0*sqr(M_PI));

    phase = (Phase*) _phase;

    if ( !(fstruct >> type >> x >> y >> z >> dw) ) return 1;

    pos[0] = x; pos[1] = y; pos[2] = z;

    // ------ - We have DISCUS format

    if (ldiscus)
    {
        u[0] = fac*dw;
        u[1] = fac*dw;
        u[2] = fac*dw;

        u[3] = u[4] = u[5] = 0.0;

        occ = 1.0;
        docc = 0.0;

        for (j=0; j<3; j++)
            dpos[j] = 0.0;

        for (j=0; j<6; j++)
            du[j] = 0.0;
    }
    else

        // ------ - We have PDFFIT format

    {
        occ = dw;
        read3(fstruct, dpos);
        fstruct >> docc;
        read3(fstruct, u);
        read3(fstruct, du);
        read3(fstruct, u+3);
        read3(fstruct, du+3);
    }

    iscat = phase->get_iscat(type);

    if (iscat == -1)
    {
        phase->at_lis.push_back(type);
        phase->Nscatlen.push_back(vector<double>());
        phase->Xscatfac.push_back(vector<double>());
        iscat = phase->nscat++;
    }

    return 0;
}

// look for atom type in the phase
int Phase::get_iscat(string type)
{
    int j;

    //do_cap(type);
    //_strupr(type);

    for (j=0; j<nscat; j++)
    {
        if (type == at_lis[j])
            return j;
    }
    return -1;
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
    const double fac = 8.0*pi*pi /3.0;
    bool ldis=false;

    // -- Write new type of structure file

    fout << "title  " << name << endl;

    fout << setprecision(6) << fixed;

    if(!ldis)
    {
        fout << "format pdffit" << endl;
        fout << "scale  " << setw(9) << skal << endl;
        fout << "sharp  " << setw(9) << delta << ", " << setw(9) << gamma << ", "
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
    for (int ia=0; ia<natoms; ia++)
    {
        Atom &atom=this->atom[ia];

        if (ldis)
        {
            double dw = fac*(atom.u[1]+atom.u[2]+atom.u[3]);
            fout << setw(4) << left << at_lis[atom.iscat];
            fout << right << setprecision(8);
            for (int i=0; i<3; i++) fout << setw(18) << atom.pos[i];
            fout << setw(13) << dw << endl;
        }
        else
        {
            fout << setw(4) << left << at_lis[atom.iscat];
            fout << right << setprecision(8);
            for (int i=0; i<3; i++) fout << setw(18) << atom.pos[i];
            fout << setw(13) << setprecision(4) << atom.occ << endl;

            fout << "    ";
            fout << setprecision(8);
            for (int i=0; i<3; i++) fout << setw(18) << atom.dpos[i];
            fout << setw(13) << setprecision(4) << atom.docc << endl;

            fout << "    ";
            fout << setprecision(8);
            for (int i=0; i<3; i++) fout << setw(18) << atom.u[i];
            fout << endl;

            fout << "    ";
            for (int i=0; i<3; i++) fout << setw(18) << atom.du[i];
            fout << endl;

            fout << "    ";
            for (int i=3; i<6; i++) fout << setw(18) << atom.u[i];
            fout << endl;

            fout << "    ";
            for (int i=3; i<6; i++) fout << setw(18) << atom.du[i];
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
       << at_lis[atom[ia].iscat] << " (#" << ia+1 << ")"
       <<" - "
       << at_lis[atom[ja].iscat] << " (#" << ja+1 << ")"
       <<" - "
       << at_lis[atom[ka].iscat] << " (#" << ka+1 << ")"
       << "   =   "  << putxdx(ang,dang) << " degrees\n";

   return ang;
}

/***************************************
c   Calculate bond lengths with errors
****************************************/
double PdfFit::bond_length(int ia, int ja)
{
    if( !curphase )
    {
        throw unassignedError("Must read structure first");
        return 0;
    }
    return curphase->bond_length(ia, ja);
}

double Phase::bond_length(int ia, int ja)
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
    cout << "   " << at_lis[atom[ia].iscat] << " (#" << ia+1 << ")" << " - "
    << at_lis[atom[ja].iscat] << " (#" << ja+1 << ")   =   " << putxdx(dist,ddist) << " A" << endl;
    return dist;
}


void PdfFit::bond_length(int ia, int ja, double bmin, double bmax)
{
    if( !curphase )
    {
        throw unassignedError("Must read structure first");
    }
    curphase->bond_length(ia, ja, bmin, bmax);
}

void Phase::bond_length(int itype, int jtype, double bmin, double bmax)
{
    bool lfound;
    double d[3], dd[3], dist, ddist;
    vector<bool> i_allowed(nscat), j_allowed(nscat);

    lfound = false;

    // ---- Get all bonds in specified range

    // get the requested atom types
    get_atoms(itype,i_allowed);
    get_atoms(jtype,j_allowed);

    cout << "(" << ((itype != ALL) ? at_lis[itype-1] : "ALL");
    cout << "," << ((jtype != ALL) ? at_lis[jtype-1] : "ALL") << ")";
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

    for (int ia=0; ia<natoms; ia++)
    {
        if (!i_allowed[ atom[ia].iscat ]) continue;
	for (int ja = 0; ja < natoms; ja++)
	{
	    if (!j_allowed[ atom[ja].iscat ])  continue;
	    for (sph.rewind(); !sph.finished(); sph.next())
	    {
		for (int jj=0; jj<3; jj++)
		{
		    d[jj] = atom[ia].pos[jj] - atom[ja].pos[jj] - 
			    sph.mno[jj]*icc[jj];
		    dd[jj] = atom[ia].dpos[jj] + atom[ja].dpos[jj];
		}
		dist = sqrt(skalpro(d,d));
		if ( (dist >= bmin) && (dist <= bmax) )
		{
		    ddist = 0.5/dist * dskalpro(d,d,dd,dd);
		    cout << "   " << at_lis[atom[ia].iscat] << " (#" << ia+1 << ")" << " - "
			<< at_lis[atom[ja].iscat] << " (#" << ja+1 << ")   =   "
			<< putxdx(dist,ddist) << " A" << endl;
		    lfound = true;
		}
	    }
	}
    }

    cout << endl;

    if (!lfound) cout << "   *** No pairs found ***\n";
}

/*
    Gets the selected atoms
*/
// ia is in [1,nscat]
void Phase::get_atoms(int ia, vector<bool> &latom)
{
    for (int i=0; i<nscat; i++)
      latom[i] = false;

    if(ia == ALL)
    {
        for (int i=0; i<nscat; i++)
            latom[i] = true;
    }
    else
    {
        if ( (ia>=1) && (ia<=nscat)) latom[ia-1] = true;
        else {
            //warning("Incorrect atom type");
            throw ValueError("Incorrect atom type");
            return;
        }
    }
}

#if FORTRAN

*****7****************************************************************
	subroutine update_cr_dim
c-
c	Updates the crystal dimensions to the current values
c+
	implicit   none
c
	include    'config.inc'
	include    'crystal.inc'
c
	integer    i,j
c
c------	Set initial values
c
	if(cr_natoms(cr_ipha).gt.0) then
	  do j=1,3
	    cr_dim(j,1,cr_ipha) = cr_pos(j,1,cr_ipha)
	    cr_dim(j,2,cr_ipha) = cr_pos(j,1,cr_ipha)
	  enddo
	else
	  do j=1,3
	    cr_dim(j,1,cr_ipha) =  1.e10
	    cr_dim(j,1,cr_ipha) = -1.e10
	  enddo
	endif
c
c------	Update values from all atoms in crystal
c
	do i=1,cr_natoms(cr_ipha)
	  do j=1,3
            cr_dim(j,1,cr_ipha)=amin1(cr_dim(j,1,cr_ipha),
     &	                              cr_pos(j,i,cr_ipha))
            cr_dim(j,2,cr_ipha)=amax1(cr_dim(j,2,cr_ipha),
     &	                              cr_pos(j,i,cr_ipha))
	  enddo
	enddo
c
c------ Set cr_dim0 values
c
	do i=1,3
	  cr_dim0(i,1,cr_ipha) = float(nint(cr_dim(i,1,cr_ipha)))
	  cr_dim0(i,2,cr_ipha) = float(nint(cr_dim(i,2,cr_ipha)))
	enddo
c
	end
c*****7****************************************************************
	subroutine do_occ (zeile,lp)
c-
c	Setting of the occupancies for atom/site
c+
	implicit      	none
c
	include 	'config.inc'
	include 	'crystal.inc'
	include 	'errlist.inc'
	include		'charact.inc'
c
	integer       	maxw
	parameter    	(maxw=20)
c
	character*(*)	zeile
	integer		lp
c
	character*200 	cpara(maxw)
	integer       	lpara(maxw),length
	integer       	ianz,iianz,is,i,j,k
	real          	werte(maxw)
	real          	wwerte(maxw)
c
	call get_params (zeile,ianz,cpara,lpara,maxw,lp)
	if (ier_num.ne.0) return
c
	if     (ianz.eq.2) then
	  iianz = ianz-1
	  call get_iscat(iianz,cpara,lpara,wwerte,maxw)
	  if(ier_num.ne.0) return
	  call del_params(1,ianz,cpara,lpara,maxw)
	  call ber_params(ianz,cpara,lpara,werte,maxw)
	  if (ier_num.ne.0) return
	  if (werte(1).ge.0.0 .and. werte(1).le.1.0) then
	    if(wwerte(1).eq.-1) then
	      do i=1,cr_natoms(cr_ipha)
	        cr_occ(i,cr_ipha) = werte(1)
	      enddo
	    else
	      do i=1,iianz
	        is = nint(wwerte(i))
 	        do k=1,cr_natoms(cr_ipha)
	          if(cr_iscat(k,cr_ipha).eq.is) then
	            cr_occ(k,cr_ipha) = werte(1)
	          endif
	        enddo
	      enddo
	    endif
	  else
	    ier_num = -103
	    ier_typ = ER_APPL
	  endif
	else
	  ier_num = -6
	  ier_typ = ER_COMM
	endif
c
	end
c*****7*****************************************************************
	subroutine do_temp (zeile,lp)
c-
c	Setting of the thermal factors for atom/site
c+
	implicit      	none
c
	include 	'config.inc'
	include 	'crystal.inc'
	include 	'errlist.inc'
	include		'charact.inc'
c
	integer       	maxw
	parameter    	(maxw=20)
c
	character*(*)	zeile
	integer		lp
c
	character*200 	cpara(maxw)
	integer       	lpara(maxw),length
	integer       	ianz,iianz,is,i,j,k
	real          	werte(maxw)
	real          	wwerte(maxw)
c
	call get_params (zeile,ianz,cpara,lpara,maxw,lp)
	if (ier_num.ne.0) return
c
	if (ianz.eq.2 .or. ianz.eq.4 .or. ianz.eq.7) then
	  iianz = ianz-1
	  call get_iscat(iianz,cpara,lpara,wwerte,maxw)
	  if(ier_num.ne.0) return
	  call del_params(1,ianz,cpara,lpara,maxw)
	  call ber_params(ianz,cpara,lpara,werte,maxw)
	  if (ier_num.ne.0) return
c
	  if (ianz.eq.1) then
	    werte(2) = werte(1)
	    werte(3) = werte(1)
	    werte(4) = 0.0
	    werte(5) = 0.0
	    werte(6) = 0.0
	  elseif (ianz.eq.3) then
	    werte(4) = 0.0
	    werte(5) = 0.0
	    werte(6) = 0.0
	  endif
c
	  if(wwerte(1).eq.-1) then
	    do i=1,cr_natoms(cr_ipha)
	      cr_u(1,i,cr_ipha) = werte(1)
	      cr_u(2,i,cr_ipha) = werte(2)
	      cr_u(3,i,cr_ipha) = werte(3)
	      cr_u(4,i,cr_ipha) = werte(4)
	      cr_u(5,i,cr_ipha) = werte(5)
	      cr_u(6,i,cr_ipha) = werte(6)
	    enddo
	  else
	    do i=1,iianz
	      is = nint(wwerte(i))
 	      do k=1,cr_natoms(cr_ipha)
	        if(cr_iscat(k,cr_ipha).eq.is) then
	          cr_u(1,k,cr_ipha) = werte(1)
	          cr_u(2,k,cr_ipha) = werte(2)
	          cr_u(3,k,cr_ipha) = werte(3)
	          cr_u(4,k,cr_ipha) = werte(4)
	          cr_u(5,k,cr_ipha) = werte(5)
	          cr_u(6,k,cr_ipha) = werte(6)
	        endif
	      enddo
	    enddo
	  endif
	else
	  ier_num = -6
	  ier_typ = ER_COMM
	endif
c
	end
c*****7*****************************************************************
c**********************************************************************
	subroutine save_atoms(ipha,strucfile,lp)
c+
c	This subroutine saves the structure and/or the unit cell
c	onto a file. Format is for plotting with ATOMS.
c+
	implicit	none
c
	include 	'config.inc'
	include 	'crystal.inc'
	include 	'errlist.inc'
c
	character*(*) 	strucfile
	integer       	ipha,lp,ist,i,j
	logical       	lread
c
	integer       	len_str
c
	data 	ist 	/17/
c
	lread = .false.
	call oeffne(ist,strucfile,'unknown',lread)
	if(ier_num.eq.0) then
c
	  write (ist, 500) cr_name(ipha)(1:len_str(cr_name(ipha)))
	  write (ist, 510) (cr_a0(i,ipha)*cr_icc(i,ipha),i=1,3),
     &	                   (cr_win(i,ipha),i=1,3)
	  write (ist, 520) 'P1'
	  write (ist, 530)
c
	  do i=1,cr_natoms(ipha)
	    write(ist,1000) cr_at_lis(cr_iscat(i,ipha),ipha),
     &	                    cr_iscat(i,ipha),
     &	                   (cr_pos(j,i,ipha)/cr_icc(j,ipha),j=1,3),
     &	                   (cr_u(j,i,ipha),j=1,6)
	  enddo
	endif
	close(ist)
c
500	format ('TITL ',a)
510	format ('CELL ',6(f10.6,1x))
520	format ('SPGP ',a)
530	format ('FIELDS LAB TYP COO TFU')
1000	format (a4,1x,i4,3x,3(f11.6,1x),/,6(f8.6,1x))
	end
c*****7*****************************************************************

#endif
