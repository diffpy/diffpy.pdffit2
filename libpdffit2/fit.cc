// Up to date with 1.3.10 Fortran version

/**********************************************************************

    This file contains all routines dealing with the
    least square refinement.

    Version  : 1.0
    Date     : 16. March 1999
    Author   : Th. Proffen (tproffen@lanl.gov)

**********************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>
#include "pdffit.h"

template <class Type> double mod(double a, Type b)
{
    int imod = int(a/b);
    if (a<0) imod--;
    return a - imod*b;
}



/*************************
    Main fit routine
**************************/
// Command <run> called <do_fit> in Fortran
// Now called <refine>

//int PdfFit::refine(bool deriv)
//{
//  double alambda, chisq, ochisq;
//  int stagnating=0;
//  const int MAXITER=100, MINITER=3, NSTAG=3;
//    //string filebase = "pdftemp";
//    //stringstream filenamestream;
//
//  // TRY OUT THE USER ROUTINES AND CATCH FOR ERRORS
//
//
//  //y: observed pdf
//  //sig: standard deviation on y (1/sqrt(wic))
//  //x: r-abscissa
//  //a: vector of parameters
//  //ia: gives possibility to fix parameters
//  //covar: returned covariance matrix
//  //alpha: returned curvature matrix
//
//  cout
//      << "*******************\n"
//      << "Starting refinement\n"
//      << "*******************\n";
//
//  for (int is=0; is<nset; is++)
//  {
//      cout << " Dataset: " << set[is]->iset << "   Phase: ";
//      for (unsigned int ip=0; ip<set[is]->psel.size(); ip++)
//          if (set[is]->psel[ip]) cout << phase[ip]->iphase << "  ";
//      cout << endl;
//
//  }
//
//
//  //------ Here starts the fitting
//
//  //_pp(fit.var); _pp(fit.ip); _pp(fit.covar.getrows()); _pp(fit.covar.getcols());
//
//  // setting the offset for all refinable variables
//
//  try
//  {
//      fit_setup();
//
//      alambda=-1;
//      fit.iter = 0;
//      do {
//          ochisq = chisq;
//
//          mrqmin(fit.p, fit.ip, fit.covar, fit.alpha, chisq, alambda, deriv);
//
//
//
//          if (fit.iter && ((ochisq-chisq) <= 1e-8*ochisq) ) stagnating++;
//          else stagnating = 0;
//          fit.iter++;
//
//
//            // save the pdf after each refinement step
//            //for( int i = 0; i < nset; i++ ) {
//            //    filenamestream << flush;
//            //    filenamestream << filebase << "_set" << i << "_iter" << fit.iter << ".gr";
//            //    save_pdf( i, filenamestream.str() );
//            //}
//
//          cout << "\n******************************** ITER: " << fit.iter << " ********************************\n";
//
//          fit.fit_rw = sqrt(chisq/fit.wnorm);
//          fit.redchisq = chisq/(fit.ntot-fit.ndof);
//
//          fit.out();
//
//          cout << " chisq.: " << chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw;
//          if (stagnating) cout << "    stagnating";
//          cout << endl;
//
//      } while( (fit.iter<MINITER) || ((stagnating < NSTAG) && (fit.iter<MAXITER)) );
//
//      cout << "\n================================ FINAL =================================\n";
//
//      alambda =0;
//      mrqmin(fit.p, fit.ip, fit.covar, fit.alpha, chisq, alambda, deriv);
//
//      fit.fit_rw = sqrt(chisq/fit.wnorm);
//      fit.redchisq = chisq/(fit.ntot-fit.ndof);
//
//      fit.out();
//
//      cout << " chisq.: " << chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw << endl;
//
//      cout << "\n======================================================================\n\n";
//
//      // final recalculation pdf with best parameters (no need for derivatives)
//      fit_theory(false,false);  // yields pdftot
//
//      fit_errors();
//
//      return 1;
//  }
//  catch(Exception e)
//  {
//      //e.PrintException();
//        throw;
//      return 0;
//  }
//}

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

    const int NSTAG = 3, MINITER = 3, MAXITER = 100;
    // Check to see if it is the first iteration.
    // If so, then set up the fit.
    // fit.iter is 0 before any refinements have been done
    // after a refinement is finished it is set to -1
    if( fit.iter <= 0 )
    {
        fit.iter = 0;
        fit.alambda = -1;
        fit.stagnating = 0;
        fit.chisq = 100;

        cout
            << "*******************\n"
            << "Starting refinement\n"
            << "*******************\n";

        for (int is=0; is<nset; is++)
        {
            cout << " Dataset: " << set[is]->iset << "   Phase: ";
            for (unsigned int ip=0; ip<set[is]->psel.size(); ip++)
                if (set[is]->psel[ip]) cout << phase[ip]->iphase << "  ";
            cout << endl;

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

        cout << "\n******************************** ITER: " << fit.iter << " ********************************\n";

        fit.fit_rw = sqrt(fit.chisq/fit.wnorm);
        fit.redchisq = fit.chisq/(fit.ntot-fit.ndof);

        fit.out();

        cout << " chisq.: " << fit.chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw;
        if (fit.stagnating) cout << "    stagnating";
        cout << endl;

        return 0;

    } else
    {


        cout << "\n================================ FINAL =================================\n";

        fit.alambda =0;
        mrqmin(fit.p, fit.ip, fit.covar, fit.alpha, fit.chisq, fit.alambda, deriv);

        fit.fit_rw = sqrt(fit.chisq/fit.wnorm);
        fit.redchisq = fit.chisq/(fit.ntot-fit.ndof);

        fit.out();

        cout << " chisq.: " << fit.chisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw << endl;

        cout << "\n======================================================================\n\n";

        // final recalculation pdf with best parameters (no need for derivatives)
        fit_theory(false,false);  // yields pdftot

        fit_errors();
        // Prepare for further fittings
        fit.iter = -1;

        return 1;
    }
}


void Fit::out()
{
    int i, j;

    cout << endl << " Refinement parameters :\n";

    for (i=0, j=0; i<psize(); i++)
    {
        if (ip[i])
        {
            cout << setw(4) << id[i] << ": " << setw(9) << fixed << p[i];

            j++;
            if (j%4) cout << "  ";
            else cout << endl;
        }
    }
    if (j%4) cout << endl;
    cout << endl;
    cout.unsetf(ios_base::fixed);
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
        fit.refvar.push_back(fit.vfind(phase.delta));
        fit.sdptr.push_back(&phase.ddelta);

        fit.refvar.push_back(fit.vfind(phase.gamma));
        fit.sdptr.push_back(&phase.dgamma);

        fit.refvar.push_back(fit.vfind(phase.skal));
        fit.sdptr.push_back(&phase.dskal);

        fit.refvar.push_back(fit.vfind(phase.srat));
        fit.sdptr.push_back(&phase.dsrat);

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
        DataSet &set=*this->set[is];

        set.offset = fit.refvar.size();

        fit.refvar.push_back(fit.vfind(set.skal));
        fit.sdptr.push_back(&set.dskal);

        fit.refvar.push_back(fit.vfind(set.sigmaq));
        fit.sdptr.push_back(&set.dsigmaq);

        fit.refvar.push_back(fit.vfind(set.qalp));
        fit.sdptr.push_back(&set.dqalp);
    }

    // maximum number of refinable variables
    int maxvar = fit.refvar.size();

    // fill the vector pointing from available refinable variables to
    // actual constrained variables

    if (maxvar != getnpar() )
    {
        //_pp(maxvar);
        //_pp(getnpar());
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

    if((this->set).size() == 0)
    {
        throw unassignedError("Space for calculation must be alloc'ed first");
        return;
    }
    //try
    //{
    fit_theory(false, true);
    return;
    //}
    //catch(Exception e)
    //{
    //  //e.PrintException();
    //  return;
    //}

}

/***********************************************************
    Here we calculate PDF and derivatives during LS fit
       (previously known as fit_theory)
************************************************************/
void PdfFit::fit_theory (bool ldiff, bool lout)
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
                phase[ip]->atom[ia].pos[i] = mod(phase[ip]->atom[ia].pos[i], phase[ip]->icc[i]);
            }
        }
    }

    for (ip=0; ip<nphase; ip++)
        phase[ip]->lattice(false);

    // determine pdf for each dataset
    fit.ntot = 0;
    fit.wnorm = 0.0;

    for (is=0; is<nset; is++)
    {
        DataSet &set=*this->set[is];

        set.determine(ldiff, lout, fit);

        //_pp(set.rmin); _pp(set.rmax); _pp(set.rfmin); _pp(set.rfmax); _pp(set.rcmin); _pp(set.rcmax);
        //_pp(set.nfmin); _pp(set.nfmax); _pp(set.ncmin); _pp(set.ncmax);

        // compute variables for reduced chi-squared and Rw
        fit.ntot += set.nfmax - set.nfmin + 1;
        for (int i=set.nfmin; i<=set.nfmax; i++)
            fit.wnorm += set.wic[i] * sqr(set.obs[i]);
    }
}

/***********************************************************************
    This routine sets up the matrix A with the derivatives d(PDF)/dx
************************************************************************/

void DataSet::fit_setup_derivatives(Fit &fit)
{
    int is, i, j, ncc, ia, ipar, offset;
    unsigned int ip;
    double fac, facs, facp, ddrho;
    double r, bk;
    //_p("Entering <fit_setup_derivatives>");
    fac = facs = facp = ddrho = r = bk = 0;

    DataSet &set=*this;

    //------ Loop over all data points

    for (i=set.ncmin; i<=set.ncmax; i++)
    {
        // --- Some common variables

        r = i*set.deltar + set.rmin;

        if (set.sigmaq > 0.0)
            bk = exp(-sqr(r*set.sigmaq)/2.0);
        else
            bk = 1.0;

        //------ ----------------------------------------------------------------
        //------     Derivatives per atom : x,y,z,u,o
        //------ ----------------------------------------------------------------

        for (ip=0; ip<psel.size(); ip++)
        {
            // only fill derivatives if phase has been selected for dataset
            if (!set.psel[ip]) continue;

            Phase &phase=*psel[ip];

            facp = phase.skal*set.skal*bk;
            facs = 1.0 / (phase.np*r);
            fac  = facs*facp;
            ncc  = phase.icc[0]*phase.icc[1]*phase.icc[2];

            for (ia=0; ia<phase.natoms; ia++)
            {
                Atom &atom=phase.atom[ia];

                offset = atom.offset;

                // -- d/dx[ip,ia], d/du[j,ip,ia]

                for (j=0;j<9; j++)
                {
                    if ( (ipar=fit.refvar[offset++]) != -1)
                        set.fit_a[i][ipar] *= fac;
                }

                // ----- ------- d/d occ[ip,ia]

                if ( (ipar=fit.refvar[offset++]) != -1)
                {
                    is = atom.iscat;

                    set.fit_a[i][ipar] = facp*(facs*set.fit_a[i][ipar]
                                    - fpi*r*phase.rho0*phase.dnorm/phase.np);

                    // all occupancies occur in np and <b>, so every contribution to
                    // the pdf contributes to the derivatives
                    set.fit_a[i][ipar] += phase.skal*set.skal*(1.0-2.0*phase.weight[is])/phase.np
                                *(calc[i][ip] + fpi*r*phase.rho0*phase.dnorm);

                }
            }

        //------ ----------------------------------------------------------------
        //-----     Derivatives per phase : lat,delta,csca,srat
        //----- ----------------------------------------------------------------

            // ----- ----- d/d(lat[j] for j=1,2,3)

            offset = phase.offset;

            for (j=0;j<3;j++)
            {
                if ( (ipar=fit.refvar[offset++]) != -1)
                    set.fit_a[i][ipar] = facp*(facs*set.fit_a[i][ipar] +
                     fpi*r*phase.dnorm*phase.rho0/phase.a0[j]);
            }

            // ----- ----- d/d(lat[4])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sina*(phase.cosa - phase.cosb*phase.cosg);
                set.fit_a[i][ipar] = facp*(facs*set.fit_a[i][ipar] +
                    fpi*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(lat[5])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sinb*(phase.cosb - phase.cosa*phase.cosg);
                set.fit_a[i][ipar] = facp*(facs*set.fit_a[i][ipar] +
                    fpi*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(lat[6])

            if ( (ipar=fit.refvar[offset++]) != -1)
            {
                ddrho = sqr(phase.a0[0]*phase.a0[1]*phase.a0[2]/phase.v)*
                    rad*phase.sing*(phase.cosg - phase.cosa*phase.cosb);
                set.fit_a[i][ipar] = facp*(facs*set.fit_a[i][ipar] +
                    fpi*r*phase.dnorm*phase.rho0*ddrho);
            }

            // ----- ----- d/d(delta[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                set.fit_a[i][ipar] *= fac;


            // ----- ----- d/d(gamma[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                set.fit_a[i][ipar] *= fac;

            //----- ----- d/d(csca[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                set.fit_a[i][ipar] = set.calc[i][ip] * set.skal;

            // ----- ----- d/d(srat[ip])

            if ( (ipar=fit.refvar[offset++]) != -1)
                set.fit_a[i][ipar] *= fac;
        }

        //------ ----------------------------------------------------------------
        //------     Derivatives per dataset : dsca,qsig
        //------ ----------------------------------------------------------------

        offset = set.offset;

        // ----- --- d/d(dsca[is])

        if ( (ipar=fit.refvar[offset++]) != -1)
            set.fit_a[i][ipar] = set.pdftot[i] / set.skal;

        // ----- --- d/d(qsig[is])

        if ( (ipar=fit.refvar[offset++]) != -1)
        {
            if (set.sigmaq > 0.0)
                set.fit_a[i][ipar] = -r*r * set.sigmaq * set.pdftot[i];
            else
                set.fit_a[i][ipar] = 0;
        }

        // ----- --- d/d(qalp[ip])

        if ( (ipar=fit.refvar[offset++]) != -1)
            set.fit_a[i][ipar] *= fac;
    }

//------ Finally we need to convolute the derivatives with
//------ the SINC function

    if (set.qmax > 0.0)
    {
        vector<double> ppp(set.ncmax+1);

        for(ip=0; ip<fit.var.size(); ip++)
        {
            if (!fit.vref[ip]) continue;

            // SHOULD THE LIMITS NOT BE NCMIN AND NCMAX AS IN DETERMINE??
            // AND WHAT WHEN NFMIN != 0 ????
            //_pp(set.ncmin); _pp(set.ncmax); _pp(set.nfmin); _pp(set.nfmax);
            // if(set.nfmin) throw "nfmin != 0";
            for(i=set.ncmin; i<=set.ncmax; i++)
            {
                int k;

                ppp[i] = set.fit_a[i][ip]*(set.qmax-set.sinc[2*i+1]);

                for (k=set.ncmin;k<=i-1;k++)
                    ppp[i] += set.fit_a[k][ip]*(set.sinc[i-k-1]-set.sinc[i+k+1]);

                for (k=i+1;k<=set.ncmax;k++)
                    ppp[i] += set.fit_a[k][ip]*(set.sinc[k-i-1]-set.sinc[k+i+1]);

            }

            for (i=set.ncmin;i<=set.ncmax;i++)
                set.fit_a[i][ip] = ppp[i];
        }
    }

    for(ip=0; ip<fit.var.size(); ip++)
    {
        if (!fit.vref[ip]) continue;

        for (i=set.ncmin;i<=set.ncmax;i++)
            set.fit_a[i][ip] *= set.deltar/pi;
    }

    // compute derivatives wrt parameters

    fit_b = fit_a * fit.dvdp;

    //_pp(set[0].fit_a[100]);
    //_pp(fit.dvdp);
    //_pp(fit_b[200]);

    //_p("Exiting <fit_setup_derivatives>");
}

#include <algorithm>

// returns the constraint index or -1 if variable not constrained
int Fit::vfind(double &a)
{
    /*vector<double*>::iterator iter;

    iter = find(var.begin(), var.end(), &a);
    int index = iter-var.begin();
    if (index < var.size()) return(index);
    else return -1;*/

    for (unsigned int i=0; i<var.size(); i++)
    {
        // when variable found, return its index in array of constraints,
        // unless it is a fixed constraint
        if(var[i] == &a)
        {
            if (vref[i])
                return i;
            else
                return -1;
        }

    }
    // if fallen through loop: variable not to be refined
    return -1;
}

void Fit::constrain(double &a, string inpform, fcon f, int ipar, FCON type)
{
    int ivar;
    //Thu Nov 17 2005 -- CLF
    //Removed this block. I don't remember when I added it. It is wrong.
    //if (!a) {
    //  throw unassignedError("Variable does not exist");
    //}
    if ( (ivar=vfind(a)) != -1)
    {
        form[ivar] = inpform;
        fconstraint[ivar] = f;
        idef[ivar] = ipar;
        ctype[ivar] = type;
        vref[ivar] = true;
        cout << "Warning: replacing existing constraint\n\n";
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
    //_p("Entering <fill_variables>");

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
                    //_pp(id[ip]); _pp(dvdp[i][ip]);
                }
            }
        }
        else if (f)
        {
            //try
            //{
            *var[i] = f(p, dvdp[i]);
            for(int ipar=0; ipar<psize(); ipar++)
                vref[i] = vref[i] || (dvdp[i][ipar] && this->ip[ipar]);
            //}
            //catch(Exception)
            //{
            //  //cout << "Error in constraint " << i << endl;
            //    throw;
            //}
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

            vref[i] = this->ip[ipar];  // constraint is fixed if parameter ipar is fixed

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
    //_p("Exiting <fill_variables>");
}

int Fit::parfind(unsigned int j)
{
    // find the position of j in parameter indices vector id
    // return -1 if not found
    vector<unsigned int>::iterator jpos;
    jpos = find(id.begin(), id.end(), j);
    if (jpos == id.end())
    {
        return -1;
    }
    return int(jpos - id.begin());
}

void Fit::setpar(unsigned int i, double val)
{
    int ipar = parfind(i);
    if (ipar != -1)
    {
        p[ipar] = val;
    }
    else
    {
        p.push_back(val);
        ip.push_back(1);    // select refinement "ON" when par gets defined
        id.push_back(i);    // store the parameter identifier
    }
}

double Fit::getpar(unsigned int n)
{
    int ipar = parfind(n);
    if (ipar < 0)
    {
        ostringstream msg;
        msg << "Parameter " << n << " does not exist";
        throw unassignedError( msg.str() );
    }
    return p[ipar];
}

void Fit::fixpar(int n)
{
    if (n==ALL)
    {
        for (int i=0; i<psize(); i++)
            ip[i] = 0;
    }
    else
    {
        int ipar = parfind(n);

        if (ipar == -1) { cout << "Warning parameter " << n << " undefined\n\n"; }
        else ip[ipar] = 0;
    }
}

void Fit::freepar(int n)
{
    if (n==ALL)
    {
        for (int i=0; i<psize(); i++)
            ip[i] = 1;
    }
    else
    {
        int ipar = parfind(n);

        if (ipar == -1) { cout << "Warning parameter " << n << " undefined\n\n"; }
        else ip[ipar] = 1;
    }
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
    //  cout << i << " " << fit.covar[i][i] << endl;

    //_pp(fit.covar);
    //_pp(fit.dp);
    //_pp(fit.dvar);

    // convert errors on refined constraints to actual pdf-parameters
    // routine very similar to fit_setup
    fit.fill_errors();

    //------ Recalculate errors for metric tensor ..

    for (ip=0; ip<nphase; ip++)
      phase[ip]->lattice(false);
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
