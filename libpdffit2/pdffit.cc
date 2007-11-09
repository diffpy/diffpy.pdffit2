/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch, Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* PdfFit and Fit methods for implementing PDFFIT1 interpreter commands.
*
* Comments:
*
* $Id$
*
***********************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "MathUtils.h"
#include "pdffit.h"

/**********************************************************
    resets the data sets and crystal structures to empty
***********************************************************/

void PdfFit::reset()
{
    int i;

    //------ Data sets

    for (i = 0; i < nset; i++)	    delete datasets[i];
    nset = 0;
    datasets.clear();

    //------ Structure

    for (i = 0; i < nphase; i++)    delete phase[i];
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
//    builtin.insert(entry("exp10",Builtin(exp10,dexp10)));
    builtin.insert(entry("log",Builtin(log,dlog)));
//    builtin.insert(entry("log10",Builtin(log10,dlog10)));
    builtin.insert(entry("sqr",Builtin(sqr,dsqr)));
    builtin.insert(entry("cube",Builtin(cube,dcube)));
    builtin.insert(entry("sqrt",Builtin(sqrt,dsqrt)));
}

/*
    Initialization routine
*/
void PdfFit::init()   // called setup in Fortran program
{
    fit.init_builtins();
}

void PdfFit::setphase(int ip)
{
    if ((ip<1) || (ip > nphase))
    {
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
    delta2.setptr(&phase.delta2);
    delta1.setptr(&phase.delta1);
    sratio.setptr(&phase.sratio);
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
        eout << "Warning: set " << is << " undefined";
        throw unassignedError(eout.str());
        return;
    }

    DataSet* pds = this->datasets[is-1];

    curset = pds;

    dscale.setptr( &(pds->skal) );
    qdamp.setptr( &(pds->qdamp) );
    qbroad.setptr( &(pds->qbroad) );
    spdiameter.setptr( &(pds->spdiameter) );
}


vector<double> DataSet::getpdf_fit()
{
        return Pdf::getpdf_fit();
}

vector<double> PdfFit::getpdf_obs()
{
    if (!curset)
    {
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
        throw unassignedError("No data loaded");
        return 0.0;
    }
    else
        return curset->rmax;
}

map<string, vector<double> > PdfFit::getPhaseFractions()
{
    if (!curset)
    {
        const char* emsg = "Dataset not defined, unknown scattering type";
        throw unassignedError(emsg);
    }
    map<string, vector<double> > rv;
    vector< pair<double,double> > atomfractions;
    vector< pair<double,double> > cellfractions;
    vector< pair<double,double> > massfractions;
    atomfractions = curset->getAtomPhaseFractions();
    cellfractions = curset->getCellPhaseFractions();
    massfractions = curset->getMassPhaseFractions();
    size_t n = atomfractions.size();
    for (size_t i = 0; i != n; ++i)
    {
        rv["atom"].push_back(atomfractions[i].first);
        rv["stdatom"].push_back(atomfractions[i].second);
        rv["cell"].push_back(cellfractions[i].first);
        rv["stdcell"].push_back(cellfractions[i].second);
        rv["mass"].push_back(massfractions[i].first);
        rv["stdmass"].push_back(massfractions[i].second);
    }
    return rv;
}

// End of file
