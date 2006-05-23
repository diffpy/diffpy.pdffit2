/****************************************************************
 PDFFIT port to C++
 translate the Proffen interpreter commands to methods of the
 pdffit class
 **************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <map>
#include "pdffit.h"
using namespace std;


string locname[]={"Application", "Command" };

string _errmsg[]={" ", "File does not exist"};
vector<string> errmsg(_errmsg,_errmsg+1);

void throw_exception(int errnum, ErrLoc errloc)
{
    errnum = abs(errnum);
    stringstream eout;
    if ((int) errmsg.size() >= errnum)
        eout << "Error: " << errmsg[errnum-1];
    else
        eout << "Error: " << errnum << " in " << locname[errloc];
    throw calculationError(eout.str());
}

void throw_exception(int errnum, string msg)
{
    errnum = abs(errnum);
    stringstream eout;
    if ( (int) errmsg.size() >= errnum)
        eout << "Error: " << errmsg[errnum-1] << " : " << msg << endl;
    else
        eout << "Error: " << errnum << " msg: " << msg << endl;
    throw calculationError(eout.str());
}

void warning(string msg)
{
    cout << "Warning: " << msg << endl;
}

/**********************************************************
    resets the data sets and crystal structures to empty
***********************************************************/

void PdfFit::reset()
{
    int i;

    //------ Data sets

    for (i=0; i<nset; i++) delete set[i];
    nset = 0;
    set.clear();

    //------ Structure

    for (i=0; i<nphase; i++) delete phase[i];
    total = 0;
    nphase = 0;
    phase.clear();

    // ------ Fit

    fit.reset();
}

void Fit::reset()
{
    // reset all data members to initial values
    alambda = chisq = ochisq = fit_rw = redchisq = wnorm = 1.0;
    stagnating =  ntot = ndof = 1;
    iter = 0;
    // clean all arrays
    p.clear();
    dp.clear();
    id.clear();
    ip.clear();
    covar.clear();
    alpha.clear();
    var.clear();
    dvar.clear();
    vref.clear();
    vcovar.clear();
    fconstraint.clear();
    form.clear();
    idef.clear();
    ctype.clear();
    used.clear();
    dvdp.clear();
    sdptr.clear();
    refvar.clear();
}

typedef pair<string,Builtin> entry;

void Fit::init_builtins()
{
    builtin.insert(entry("-",Builtin(neg,dneg)));
    builtin.insert(entry("sin",Builtin(sin,dsin)));
    builtin.insert(entry("cos",Builtin(cos,dcos)));
    builtin.insert(entry("tan",Builtin(tan,dtan)));
    builtin.insert(entry("sind",Builtin(sind,dsind)));
    builtin.insert(entry("cosd",Builtin(cosd,dcosd)));
    builtin.insert(entry("tand",Builtin(tand,dtand)));
    builtin.insert(entry("asin",Builtin(asin,dasin)));
    builtin.insert(entry("acos",Builtin(acos,dacos)));
    builtin.insert(entry("atan",Builtin(atan,datan)));
    builtin.insert(entry("asind",Builtin(asind,dasind)));
    builtin.insert(entry("acosd",Builtin(acosd,dacosd)));
    builtin.insert(entry("atand",Builtin(atand,datand)));
    builtin.insert(entry("exp",Builtin(exp,dexp)));
    builtin.insert(entry("exp10",Builtin(exp10,dexp10)));
    builtin.insert(entry("log",Builtin(log,dlog)));
    builtin.insert(entry("log10",Builtin(log10,dlog10)));
    builtin.insert(entry("sqr",Builtin(sqr,dsqr)));
    builtin.insert(entry("cube",Builtin(cube,dcube)));
    builtin.insert(entry("sqrt",Builtin(sqrt,dsqrt)));
}

/*
    Initialization routine
*/
void PdfFit::init()   // called setup in Fortran program
{
//
//------ Write starting screen
//
    string sCreatedDate("Created : ");
    sCreatedDate += cdate;
    cout << endl;
    cout << "          ***********************************************************\n";
    cout << "          *               P D F F I T   Version       "
                          << left << setw(14) << version << "*\n";
    cout << "          *                                                         *\n";
    cout << "          *         " << right << setw(45) << sCreatedDate <<   "   *\n";
    cout << "          *---------------------------------------------------------*\n";
    cout << "          * (c) Thomas Proffen  -     Email: tproffen@lanl.gov      *\n";
    cout << "          *     Simon Billinge  -     Email: billinge@pa.msu.edu    *\n";
    cout << "          *     Jacques Bloch   -     Email: bloch@pa.msu.edu       *\n";
    cout << "          ***********************************************************\n";
    cout << endl;

    fit.init_builtins();

//
//------    get envirmonment information
//
    // appl_env();  will read environment variables: in appl_win.f and appl_unix.f
//
//------    try to read default file
//
    // autodef();
//
//------    try to read command line arguments
//
    // cmdline_args();

}

void PdfFit::setphase(int ip)
{
    if ((ip<1) || (ip > nphase))
    {
        //cout << "Warning: phase " << ip << " undefined\n";
        stringstream eout;
        eout << "Warning: phase " << ip << " undefined";
        throw unassignedError(eout.str());
        return;
    }

    Phase &phase=*this->phase[ip-1];

    curphase = &phase;

    lat.resize(6);
    lat[0].setptr(&phase.a0[0]);
    lat[1].setptr(&phase.a0[1]);
    lat[2].setptr(&phase.a0[2]);
    lat[3].setptr(&phase.win[0]);
    lat[4].setptr(&phase.win[1]);
    lat[5].setptr(&phase.win[2]);

    pscale.setptr(&phase.skal);
    delta.setptr(&phase.delta);
    gamma.setptr(&phase.gamma);
    srat.setptr(&phase.srat);
    rcut.setptr(&phase.rcut);

    x.resize(phase.natoms);
    y.resize(phase.natoms);
    z.resize(phase.natoms);
    u11.resize(phase.natoms);
    u22.resize(phase.natoms);
    u33.resize(phase.natoms);
    u12.resize(phase.natoms);
    u13.resize(phase.natoms);
    u23.resize(phase.natoms);
    occ.resize(phase.natoms);

    for(int ia=0; ia<phase.natoms; ia++)
    {
        Atom &atom=phase.atom[ia];

        x[ia].setptr(&atom.pos[0]);
        y[ia].setptr(&atom.pos[1]);
        z[ia].setptr(&atom.pos[2]);

        u11[ia].setptr(&atom.u[0]);
        u22[ia].setptr(&atom.u[1]);
        u33[ia].setptr(&atom.u[2]);
        u12[ia].setptr(&atom.u[3]);
        u13[ia].setptr(&atom.u[4]);
        u23[ia].setptr(&atom.u[5]);

        occ[ia].setptr(&atom.occ);
    }
}

void PdfFit::setdata(int is)
{
    if ((is<1) || (is > nset))
    {
        stringstream eout;
        //cout << "Warning: set " << is << " undefined\n";
        eout << "Warning: set " << is << " undefined";
        throw unassignedError(eout.str());
        return;
    }

    DataSet &set=*this->set[is-1];

    curset = &set;

    dscale.setptr(&set.skal);
    qsig.setptr(&set.sigmaq);
    qalp.setptr(&set.qalp);
}


vector<double> DataSet::getpdf_fit()
{
        return Pdf::getpdf_fit();
}

vector<double> PdfFit::getpdf_obs()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        vector <double> empty_double;
        return empty_double;
    }
    else
        return curset->getpdf_obs();
}

vector<double> PdfFit::getpdf_fit()
{
    if (!curset)
    {
        //warning("No fit yet!");
        throw unassignedError("No fit data");
        vector<double> empty_double;
        return empty_double;
    }
    else
        return curset->getpdf_fit();
}

int PdfFit::getnfmin()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        return 0;
    }
    else
        return curset->nfmin;
}

int PdfFit::getnfmax()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        return 0;
    }
    else
        return curset->nfmax;
}

double PdfFit::getdeltar()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        return 0.0;
    }
    else
        return curset->deltar;
}
double PdfFit::getrmin()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        return 0.0;
    }
    else
        return curset->rmin;
}
double PdfFit::getrmax()
{
    if (!curset)
    {
        //warning("No data loaded!");
        throw unassignedError("No data loaded");
        return 0.0;
    }
    else
        return curset->rmax;
}
