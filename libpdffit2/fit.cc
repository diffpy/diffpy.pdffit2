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
* Mixed definitions of several DataSet, Fit and PdfFit methods
*
* Comments: Up to date with 1.3.10 Fortran version.
*	    What a spagetti.
*
***********************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <valarray>

#include "MathUtils.h"
#include "ShapeFactors.h"
#include "pdffit.h"

using NS_PDFFIT2::pout;

/*************************
    Main fit routine
**************************/
// Command <run> called <do_fit> in Fortran
// Now called <refine>

int PdfFit::refine(bool deriv, double toler) {
    int finished = 0;
    while( !finished )
    {
        finished = refine_step(deriv, toler);
    }
    return 1;
}


// CLF Wed May 25 2005
// Added one-step refinement so that the progress can be checked after each
// refinement step.  Many parameters made into public data members of Fit. This
// function is for use in higher level refinement routines.
int PdfFit::refine_step(bool deriv, double toler)
{
    static bool fit_running = false;
    const int NSTAG = 3, MINITER = 3, MAXITER = 100;

    fit_running = fit_running && (fit.iter != 0);
    // If fit_running flag is down, this is the first iteration.
    // If so, then set up the fit.
    if(!fit_running)
    {
        fit.iter = 0;
        fit.alambda = -1;
        fit.stagnating = 0;
        fit.chisq = double_max;
	fit_running = true;

	*pout <<
	    "*******************\n" <<
	    "Starting refinement\n" <<
	    "*******************\n";

        for (int is=0; is<nset; is++)
        {
            *pout << " Dataset: " << datasets[is]->iset << "   Phase: ";
            for (unsigned int ip=0; ip<datasets[is]->psel.size(); ip++)
                if (datasets[is]->psel[ip]) *pout << phase[ip]->iphase << "  ";
            *pout << endl;

        }

        fit_setup();
    }


    //------ Here starts the fitting

    // setting the offset for all refinable variables

    if( (fit.iter<MINITER) || ((fit.stagnating < NSTAG) && (fit.iter<MAXITER)) )
    {
        fit.ochisq = fit.chisq;

        mrqmin(fit.p, fit.ip, fit.covar, fit.alpha, fit.chisq, fit.alambda, deriv);

        if (fit.iter && ((fit.ochisq-fit.chisq) <= toler*fit.ochisq) )
        {
            fit.stagnating++;
        }
        else fit.stagnating = 0;
        fit.iter++;

        *pout << "\n******************************** ITER: " << fit.iter << " ********************************\n";

        fit.fit_rw = sqrt(fit.chisq/fit.wnorm);
        fit.redchisq = fit.chisq/(fit.ntot-fit.ndof);

        fit.out();

        *pout << " chisq.: " << fit.chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw;
        if (fit.stagnating) *pout << "    stagnating";
        *pout << endl;

        return 0;

    }
    else
    {

        *pout << "\n================================ FINAL ==================================\n";

        fit.alambda =0;
        mrqmin(fit.p, fit.ip, fit.covar, fit.alpha, fit.chisq, fit.alambda, deriv);

        fit.fit_rw = sqrt(fit.chisq/fit.wnorm);
        fit.redchisq = fit.chisq/(fit.ntot-fit.ndof);

        fit.out();

        *pout << " chisq.: " << fit.chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw << endl;

        *pout << "\n=========================================================================\n\n";

        // final recalculation pdf with best parameters (no need for derivatives)
        fit_theory(false,false);  // yields pdftot

        fit_errors();

        // Prepare for further fittings
        fit_running = false;

        return 1;
    }
}


void Fit::out()
{
    *pout << endl << " Refinement parameters :\n";

    vector<int> order = this->order_by_id();
    vector<int>::iterator i;
    int j;
    for (i = order.begin(), j = 0; i != order.end(); ++i)
    {
        if (ip[*i])
        {
            *pout << setw(4) << id[*i] << ": " << setw(9) << fixed << p[*i];
            j++;
            if (j % 4) *pout << "  ";
            else *pout << endl;
        }
    }
    if (j % 4) *pout << endl;
    *pout << endl;
    (*pout).unsetf(ios_base::fixed);
}

/**********************************************************
    Setting the offset for all refinable variables
**********************************************************/
void PdfFit::fit_setup()
{
    // initialize the array with the addresses of refinable variable
    // and initialize the offset pointers to be used in derivative routines
    // also makes a vector containing the address of the error of each
    // refinable variable
    //int ip, ia, is, j;
    int j;

    // First we make an initial fill of all constrained equations
    // to detect missing parameters but mainly to detect the fixed constraints,
    // i.e. constraints without parameters or with only fixed parameters

    fit.fill_variables();

    // for each refinable variable detect the corresponding constraint

    fit.refvar.clear();
    fit.sdptr.clear();

    for (int ip=0; ip<nphase; ip++)
    {
        Phase &phase=*this->phase[ip];

        phase.offset = fit.refvar.size();

        for (j=0;j<3;j++)
        {
            fit.refvar.push_back(fit.vfind(phase.a0[j]));
            fit.sdptr.push_back(&phase.da0[j]);
        }
        for (j=0;j<3;j++)
        {
            fit.refvar.push_back(fit.vfind(phase.win[j]));
            fit.sdptr.push_back(&phase.dwin[j]);

        }
        fit.refvar.push_back(fit.vfind(phase.delta2));
        fit.sdptr.push_back(&phase.ddelta2);

        fit.refvar.push_back(fit.vfind(phase.delta1));
        fit.sdptr.push_back(&phase.ddelta1);

        fit.refvar.push_back(fit.vfind(phase.pscale));
        fit.sdptr.push_back(&phase.dpscale);

        fit.refvar.push_back(fit.vfind(phase.spdiameter));
        fit.sdptr.push_back(&phase.dspdiameter);

        fit.refvar.push_back(fit.vfind(phase.sratio));
        fit.sdptr.push_back(&phase.dsratio);

        for (int ia=0; ia<phase.natoms; ++ia)
        {
            Atom &atom=phase.atom[ia];

            atom.offset = fit.refvar.size();

            for (j=0;j<3; j++)
            {
                fit.refvar.push_back(fit.vfind(atom.pos[j]));
                fit.sdptr.push_back(&atom.dpos[j]);

            }
            for (j=0;j<6; j++)
            {
                fit.refvar.push_back(fit.vfind(atom.u[j]));
                fit.sdptr.push_back(&atom.du[j]);

            }
            fit.refvar.push_back(fit.vfind(atom.occ));
            fit.sdptr.push_back(&atom.docc);
        }
    }

    for (int is=0; is<nset; is++)
    {
        DataSet& ds = *this->datasets[is];

        ds.offset = fit.refvar.size();

        fit.refvar.push_back(fit.vfind(ds.dscale));
        fit.sdptr.push_back(&ds.ddscale);

        fit.refvar.push_back(fit.vfind(ds.qdamp));
        fit.sdptr.push_back(&ds.dqdamp);

        fit.refvar.push_back(fit.vfind(ds.qbroad));
        fit.sdptr.push_back(&ds.dqbroad);
    }

    // maximum number of refinable variables
    int maxvar = fit.refvar.size();

    // fill the vector pointing from available refinable variables to
    // actual constrained variables

    if (maxvar != getnpar() )
    {
        throw constraintError("Parameter set but not constrained.");
    }

    fit.ndof = 0;

    for (int i=0; i<fit.psize(); i++)
        if (fit.ip[i]) fit.ndof++;

    fit.dp.resize(fit.p.size());
    fit.covar.resize(fit.p.size(),fit.p.size());
    fit.alpha.resize(fit.p.size(),fit.p.size());

}

// user-callable routine to compute pdf from a model
void PdfFit::calc()
{
    if (this->datasets.empty())
    {
        throw unassignedError("Space for calculation must be alloc'ed first");
    }
    fit_theory(false, true);
    *pout << "\n================================== DONE =================================\n";
    return;
}

/***********************************************************
    Here we calculate PDF and derivatives during LS fit
       (previously known as fit_theory)
************************************************************/
void PdfFit::fit_theory(bool ldiff, bool lout)
{
    int is, ip, ia, i;


// ------ First we compute all constrained equations and the PDF

    fit.fill_variables();

    // reposition atoms in the elementary unit cells
    // Check if this is the correct wayto do things
    for (ip=0; ip<nphase; ip++)
    {
        for (ia=0; ia<phase[ip]->natoms; ia++)
        {
            for(i=0; i<3; i++)
            {
                phase[ip]->atom[ia].pos[i] = fmod(phase[ip]->atom[ia].pos[i], phase[ip]->icc[i]);
            }
        }
    }

    for (ip=0; ip<nphase; ip++)
        phase[ip]->lattice();

    // determine pdf for each dataset
    fit.ntot = 0;
    fit.wnorm = 0.0;

    for (is=0; is<nset; is++)
    {
        DataSet& ds = *this->datasets[is];

        ds.determine(ldiff, lout, fit);

        // compute variables for reduced chi-squared and Rw
        fit.ntot += ds.nfmax - ds.nfmin + 1;
        fit.wnorm += ds.weighedSquareObs();
    }
}


double PdfFit::totalWeighedSquareObs() const
{
    return fit.wnorm;
}


int PdfFit::totalReducedObservations() const
{
    return fit.ntot - fit.ndof;
}


/***********************************************************************
    This routine sets up the matrix A with the derivatives d(PDF)/dx
************************************************************************/

void DataSet::fit_setup_derivatives(Fit &fit)
{
    int i, j, ia, ipar, offset;
    unsigned int ip;
    double fac, facs, facp, ddrho;
    double r, bk;
    fac = facs = facp = ddrho = r = bk = 0;

    DataSet& ds = *this;

    //------ Loop over all data points

    for (i = ds.ncmin; i <= ds.ncmax; i++)
    {
        // --- Some common variables

        r = i*ds.deltar + ds.rmin;

	// background envelope due to Q resolution
	bk = (ds.qdamp > 0.0) ? exp(-sqr(r*ds.qdamp)/2.0) : 1.0;


        //------ ----------------------------------------------------------------
        //------     Derivatives per atom : x,y,z,u,o
        //------ ----------------------------------------------------------------

        for (ip=0; ip<psel.size(); ip++)
        {
            // only fill derivatives if phase has been selected for dataset
            if (!ds.psel[ip]) continue;

            Phase& phase = *psel[ip];

            double shape_env = (phase.spdiameter <= 0.0) ? (1.0) :
                sphereEnvelope(r, phase.spdiameter);

            if (phase.stepcut > 0.0 && r > phase.stepcut)
            {
                shape_env = 0.0;
            }

            facp = phase.pscale * shape_env * ds.dscale * bk;
            facs = 1.0 / (phase.np*r);
            fac  = facs*facp;

            for (ia=0; ia<phase.natoms; ia++)
            {
                Atom &atom=phase.atom[ia];

                offset = atom.offset;

                // -- d/dx[ip,ia], d/du[j,ip,ia]

                for (j=0;j<9; j++)
                {
                    if ( (ipar=fit.refvar[offset++]) != -1)
                        ds.fit_a[i][ipar] *= fac;
                }

                // ----- ------- d/d occ[ip,ia]

                if ( (ipar=fit.refvar[offset++]) != -1)
                {
                    ds.fit_a[i][ipar] = facp*(facs*ds.fit_a[i][ipar] -
                            4.0*M_PI*r*phase.rho0*phase.dnorm/phase.np);
                    // all occupancies occur in np and <b>, so every contribution to
                    // the pdf contributes to the derivatives
                    ds.fit_a[i][ipar] += phase.pscale * ds.dscale *
                        (1.0 - 2.0 * atom.weight) / phase.np *
                        (calc[i][ip] +
                         4.0*M_PI * r * phase.rho0 * phase.dnorm * bk * shape_env);
                }
            }

        //------ ----------------------------------------------------------------
        //------    Derivatives per phase : lat, delta2, pscale, spdiameter, sratio
        //------ ----------------------------------------------------------------

            // ----- ----- d/d(lat[j] for j=1,2,3)

            offset = phase.offset;

            for (j=0;j<3;j++)
            {
                if ( (ipar=fit.refvar[offset++]) != -1)
                    ds.fit_a[i][ipar] = facp*(facs*ds.fit_a[i][ipar] +
                     4.0*M_PI*r*phase.dnorm*phase.rho0/phase.a0[j]);
            }

            // ----- ----- d/d(lat[4])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sina*(phase.cosa - phase.cosb*phase.cosg);
                ds.fit_a[i][ipar] = facp*(facs*ds.fit_a[i][ipar] +
                    4.0*M_PI*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(lat[5])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sinb*(phase.cosb - phase.cosa*phase.cosg);
                ds.fit_a[i][ipar] = facp*(facs*ds.fit_a[i][ipar] +
                    4.0*M_PI*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(lat[6])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sing*(phase.cosg - phase.cosa*phase.cosb);
                ds.fit_a[i][ipar] = facp*(facs*ds.fit_a[i][ipar] +
                    4.0*M_PI*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(delta2[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                ds.fit_a[i][ipar] *= fac;


            // ----- ----- d/d(delta1[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                ds.fit_a[i][ipar] *= fac;

            //----- ----- d/d(pscale[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                ds.fit_a[i][ipar] = ds.calc[i][ip] * ds.dscale;

            // ----- --- d/d(spdiameter)
            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ds.fit_a[i][ipar] = (phase.spdiameter <= 0.0) ? 0.0 :
                    ds.calc[i][ip] * ds.dscale * phase.pscale *
                    dsphereEnvelope(r, phase.spdiameter) /
                    ((shape_env > 0.0) ? shape_env : 1.0);
            }

            // ----- ----- d/d(sratio[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                ds.fit_a[i][ipar] *= fac;
        }

        //------ ----------------------------------------------------------------
        //------     Derivatives per dataset : dscale, qdamp, qbroad
        //------ ----------------------------------------------------------------

        offset = ds.offset;

        // ----- --- d/d(dscale[is])

        if ( (ipar=fit.refvar[offset++]) != -1)
            ds.fit_a[i][ipar] = ds.pdftot[i] / ds.dscale;

        // ----- --- d/d(qdamp[is])

        if ( (ipar=fit.refvar[offset++]) != -1)
        {
            if (ds.qdamp > 0.0)
                ds.fit_a[i][ipar] = -r*r * ds.qdamp * ds.pdftot[i];
            else
                ds.fit_a[i][ipar] = 0;
        }

        // ----- --- d/d(qbroad[ip])

        if ( (ipar=fit.refvar[offset++]) != -1)
            ds.fit_a[i][ipar] *= fac;

    }

//------ Finally we need to apply Qmax cutoff on the derivatives

    if (ds.qmax > 0.0)
    {
	int nclen = ds.ncmax + 1 - ds.ncmin;
	// matrix column is not a continuous data block, a copy is required
	valarray<double> col_ip_array(nclen);
	double* col_ip = &(col_ip_array[0]);
        for(ip=0; ip<fit.var.size(); ip++)
        {
            if (!fit.vref[ip])	    continue;
	    for (i = 0; i < nclen; ++i)    col_ip[i] = ds.fit_a[ncmin+i][ip];
	    applyQmaxCutoff(col_ip, nclen);
	    for (i = 0; i < nclen; ++i)    ds.fit_a[ncmin+i][ip] = col_ip[i];
        }
    }

    // compute derivatives wrt parameters

    fit_b = fit_a * fit.dvdp;

}


// returns the constraint index or -1 if variable not constrained
int Fit::vfind(double &a)
{
    vector<double*>::iterator apos;
    apos = find(var.begin(), var.end(), &a);
    if (apos == var.end())  return -1;
    // variable is found here, now check if it is refinable
    int idx = apos - var.begin();
    return vref[idx] ? idx : -1;
}


void Fit::constrain(double &a, string inpform, fcon f, int ipar, FCON type)
{
    int ivar;
    if ( (ivar=vfind(a)) != -1)
    {
        form[ivar] = inpform;
        fconstraint[ivar] = f;
        idef[ivar] = ipar;
        ctype[ivar] = type;
        vref[ivar] = true;
        *pout << "Warning: replacing existing constraint\n\n";
    }
    else
    {
        var.push_back(&a);
        form.push_back(inpform);
        fconstraint.push_back(f);
        idef.push_back(ipar);
        ctype.push_back(type);
        vref.push_back(true);
    }
}


void Fit::constrain(double &a, string inpform)
{
    constrain(a, inpform, NULL, -1, USER);
}

void Fit::constrain(double &a, fcon f )
{
    constrain(a, string(), f, -1, USER);
}

// if a # is passed instead of a function, then the default function will be used
// which is just var = p[ipar]
void Fit::constrain(double &a, int ipar)
{
    constrain(a, ipar, IDENT);
}

// if a # and the FCOMP-type are passed instead of a function, then the
// complement function will be used:
// which is just var = 1-p[ipar]
void Fit::constrain(double &a, int ipar, FCON type)
{
    if ( (type == IDENT) || (type == FCOMP) || (type == FSQR) )
        constrain(a, string(), NULL, ipar, type);
    else
        throw constraintError("Unknown constraint");
}

void Fit::fill_variables()
{
    dvdp.resize(var.size(),p.size());

    for (unsigned int i=0; i<var.size(); i++)
    {
        fcon f = fconstraint[i];

        vref[i] = false;

        // test if formula exist for this constraint
        if (!form[i].empty())
        {
            vector<double> dnumdp;

            *var[i] = parse(form[i],dnumdp);

            // store numerical derivatives
            for(unsigned int iu=0; iu<used.size(); iu++)
            {
                int ip=used[iu];

                vref[i] = vref[i] || this->ip[ip];

                // calculate derivative wrt to p[ip] if parameter is free
                if (this->ip[ip])
                {
                    dvdp[i][ip] = dnumdp[iu];
                }
            }
        }
        else if (f)
        {
            //try
            //{
	    vector<double> dvdp_i(dvdp.rowVector(i));
            *var[i] = f(p, dvdp_i);
	    copy(dvdp_i.begin(), dvdp_i.end(), dvdp[i]);
            for(int ipar=0; ipar<psize(); ipar++)
	    {
                vref[i] = vref[i] || (dvdp[i][ipar] && this->ip[ipar]);
            }
        }
        else
        {
            int ipar = parfind(idef[i]);

            if (ipar == -1)
            {
                ostringstream msg;
                msg << "parameter " << idef[i] << " undefined";
                throw constraintError( msg.str() );
            }

            // constraint is fixed if parameter ipar is fixed
            vref[i] = this->ip[ipar];

            if (ctype[i] == IDENT)
            {
                *var[i] = p[ipar];
                dvdp[i][ipar] = 1.0;
            }
            else if (ctype[i] == FCOMP)
            {
                *var[i] = 1.0 - p[ipar];
                dvdp[i][ipar] = -1.0;
            }
            else if (ctype[i] == FSQR)
            {
                *var[i] = sqr(p[ipar]);
                dvdp[i][ipar] = 2.0*p[ipar];
            }
        }
    }
}

int Fit::parfind(unsigned int pidx)
{
    // find the position of pidx in parameter indices vector id
    // return -1 if not found
    vector<unsigned int>::iterator pos;
    pos = find(id.begin(), id.end(), pidx);
    if (pos == id.end())
    {
        return -1;
    }
    return int(pos - id.begin());
}

void Fit::setpar(unsigned int pidx, double val)
{
    int ipar = parfind(pidx);
    if (ipar != -1)
    {
        p[ipar] = val;
    }
    else
    {
        p.push_back(val);
        ip.push_back(1);    // select refinement "ON" when par gets defined
        id.push_back(pidx);    // store the parameter identifier
    }
}

double Fit::getpar(unsigned int pidx)
{
    int ipar = parfind(pidx);
    if (ipar < 0)
    {
        ostringstream msg;
        msg << "Parameter " << pidx << " does not exist";
        throw unassignedError( msg.str() );
    }
    return p[ipar];
}

void Fit::fixpar(int pidx)
{
    if (pidx == ALL)
    {
	fill(ip.begin(), ip.end(), false);
	return;
    }
    int ipar = parfind(pidx);
    if (ipar == -1)
    {
	ostringstream emsg;
	emsg << "Parameter " << pidx << " not defined.";
	throw unassignedError(emsg.str());
    }
    ip[ipar] = false;
}

void Fit::freepar(int pidx)
{
    if (pidx == ALL)
    {
	fill(ip.begin(), ip.end(), true);
	return;
    }
    int ipar = parfind(pidx);
    if (ipar == -1)
    {
	ostringstream emsg;
	emsg << "Parameter " << pidx << " not defined.";
	throw unassignedError(emsg.str());
    }
    ip[ipar] = true;
}

/**************************************************************************
    This routine calculates errors dx from dp ...

    The covariance of the constrained variables is given by the left and right
    matrix multiplication of the dvdp matrix with the covariance matrix
    (based on a Taylor expansion of each constraint around the mean parameter values).
    The old Fortran program had a sqrt(sum of squares) as error. This
    does not take into account parameter correlations, and implicitly considers
    the parameters uncorrelated, which is an under-estimation of the parameters.
****************************************************************************/
void  PdfFit::fit_errors()
{
    int ip;
    matrix<double> dvdpt = fit.dvdp.transposed();

    // compute the errors on the parameters as sqrt of diagonal elements of
    // covariance matrix
    fit.dp = fit.covar.sd();

    // compute errors on refined variables: left and right matrix multiplication
    // of covariance matrix with dvdp
    fit.vcovar = fit.dvdp * fit.covar * dvdpt;

    fit.dvar = fit.vcovar.sd();

    //for(int i=0;i<fit.psize();i++)
    //  *pout << i << " " << fit.covar[i][i] << endl;

    // convert errors on refined constraints to actual pdf-parameters
    // routine very similar to fit_setup
    fit.fill_errors();

    //------ Recalculate errors for metric tensor ..

    for (ip=0; ip<nphase; ip++)
      phase[ip]->lattice();
}

/***********************************************************************
    This routine converts the errors <dvar> on the constraints to
    errors on the actual pdf-paramters
************************************************************************/

void Fit::fill_errors()
{
    // transfers the errors on the constrained variables to the corresponding
    // pdf error-variables

    int icon;

    // loop over all refinable variables and transfer the sd on any constrained
    // variable into the sd on the refinable variable
    for (unsigned int i=0; i<sdptr.size(); i++)
    {
        if ( (icon=refvar[i]) != -1)
            *sdptr[i] = dvar[icon];
    }
}

// End of file
