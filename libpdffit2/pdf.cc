// Up to date with 1.3.10 Fortran version

#include <fstream>
#include <sstream>
#include "pdffit.h"
#include "matrix.h"
#include <iomanip>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "PointsInSphere.h"

#define NEW_SHARP

/*****************************************************************
  Allocating space for dummy dataset when calculating PDF
  without data
******************************************************************/
void PdfFit::alloc(Sctp t, double qmax, double sigmaq, double rmin, double rmax, int bin)
{

    DataSet *_set = new DataSet, &set=*_set;
    int i;
    if (rmax < rmin || rmin < 0 || rmax < 0)
    {
        throw ValueError("Check rmin, rmax");
        return;
    }
    if (bin <= 1)
    {
        throw ValueError("bin must be > 1");
        return;
    }
    if( qmax < 0.0 )
    {
        throw ValueError("qmax must be >= 0");
        return;
    }
    if( sigmaq < 0 )
    {
        throw ValueError("sigmaq must be >= 0");
        return;
    }

    //check to see if a structure has been loaded
    if((this->phase).size() == 0)
    {
        throw unassignedError("Structure must be read first");
        return;
    }


    set.iset = nset+1;

    if (t == N)
    {
        set.lxray = false;
    }
    else if (t == X)
    {
        set.lxray = true;
    }

    set.qmax   = qmax;
    set.sigmaq = sigmaq;
    set.rmin   = set.rfmin = rmin;
    set.rmax   = set.rfmax = rmax;
    set.bin    = bin;
    set.deltar = (rmax-rmin)/double(bin-1);
    set.name   = "Dummy set";

    set.obs.resize(bin);
    set.wic.resize(bin);

    for (i=0; i<bin; i++)
    {
        set.obs[i] = 0.0;
        set.wic[i] = 1.0;
    }

    cout << " Allocated PDF data set " << set.iset << "  (r = "
        << rmin << " to " << rmax << " A, " << bin << " points) ..." << endl;

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        set.selphase(ip, this->phase[ip]);

    this->set.push_back(&set);
    nset++;
    setdata(nset);
}

// cut-away high Q harmonics using fast Fourier transformation
void DataSet::applyQmaxCutoff(double* y, size_t len)
{
    // pad y with the same number of zeros up to the next power of 2
    size_t padlen = 2*len;
    padlen = size_t( pow(2, ceil(log2(padlen))) );
    // ycpad is complex, so it needs to be twice as long
    double ycpad[2*padlen];
    fill_n(ycpad, 2*padlen, 0.0);
    // copy y to real components of ycpad
    for (size_t i = 0; i != len; ++i)	ycpad[2*i] = y[i];
    // apply fft
    int status;
    status = gsl_fft_complex_radix2_forward(ycpad, 1, padlen);
    if (status != GSL_SUCCESS)
    {
        throw calculationError("Fourier Transformation failed.");
    }
    // Q step for ycpad
    double dQ = 2*M_PI/((padlen-1)*deltar);
    // loQidx, hiQidx correspond to Qmax and -Qmax frequencies
    // they need to be integer to catch cases with huge qmax/dQ
    int loQidx = int( ceil(qmax/dQ) );
    int hiQidx = padlen + 1 - loQidx;
    // zero high Q components in ycpad
    for (int i = loQidx; i < hiQidx; ++i)
    {
	ycpad[2*i] = ycpad[2*i+1] = 0.0;
    }
    // transform back
    status = gsl_fft_complex_radix2_inverse(ycpad, 1, padlen);
    if (status != GSL_SUCCESS)
    {
        throw calculationError("Fourier Transformation failed.");
    }
    // copy real components
    for (size_t i = 0; i != len; ++i)	y[i] = ycpad[2*i];
}


/*****************************************************
    Calculate PDF for dataset with current structure

      In initial version: one phase
******************************************************/

void DataSet::determine(bool ldiff, bool lout, Fit &fit)
{
    int i, ia, is, j, iatom, jj, js;
    //int kk, ibin, ip;
    int igaus, ib, ie, ig;
    double rmax2, rmin2, gaus, rk, rb,re, rtot, ract;
    vector<double> ppp;
    double dd[3], d[3], dist2, dist, rg, r;
    double sigmap, sigma, gnorm, ampl;
    bool ldone;
    long totcalc=0;

    //_p("Entering <determine>");

    ldone = false;

    if (lout) cout << " Calculating PDF ...\n";

    // total # points to be fitted
    nfmin = nint((rfmin-rmin)/deltar);
    nfmax = nint((rfmax-rmin)/deltar);

    // extend calculation range before applying Qmax cutoff
    extendCalculationRange(lout);

    rmax2 = sqr(rcmax+1.0);
    rmin2 = sqr(rcmin-1.0) * ((rcmin-1.0) < 0.0 ? -1.0 : 1.0);
    rmin2 = max(0.01,rmin2);

    pdftot.resize(ncmax+1);
    fill(pdftot.begin(), pdftot.end(), 0.0);
    if (ncmax + 1 > int(pdfref.size()))	    pdfref.resize(ncmax+1, 0.0);
    calc.resize(ncmax+1,psel.size());
    ppp.resize(ncmax+1);

    if (ldiff)
    {
        fit_a.resize(ncmax+1, fit.varsize());
	fit_a = 0.0;
    }

    // ------ Loop over all phases

    for (unsigned ip=0; ip<psel.size(); ip++)
    {
        if (!psel[ip]) continue;  // skip not-selected phases

        Phase &phase=*psel[ip];

        phase.setup_weights(false,lxray);


        // ------ Get the ratio total pairs/selected weight in structure

        rtot = 0.0;
        ract = 0.0;

        for (i=0; i<phase.natoms; i++)
        {
            for (j=0; j<phase.natoms; j++)
            {
                is = phase.atom[i].iscat;
                js = phase.atom[j].iscat;
                rtot += phase.weight[is]*phase.weight[js];
                if ( (allowed_i[ip][is] && allowed_j[ip][js]) || (allowed_j[ip][is] && allowed_i[ip][js]) )
                    ract += phase.weight[is]*phase.weight[js];
            }
        }

        phase.dnorm = ract/rtot;

        //--------------------------------------------------------------

        for (i=ncmin; i<=ncmax; i++)
            ppp[i] = 0.0;

	// calculate range for PointsInSphere sequencer
	// (negative rsphmin is no problem)
	double buffzone = phase.circum_diameter();
	double rsphmin = sqrt(rmin2) - buffzone;
	double rsphmax = sqrt(rmax2) + buffzone;
	PointsInSphere sph( rsphmin, rsphmax, phase.a0[0]*phase.icc[0],
		phase.a0[1]*phase.icc[1], phase.a0[2]*phase.icc[2],
		phase.win[0], phase.win[1], phase.win[2] );
        for (ia=0; ia<phase.natoms; ia++)
        {
            Atom atomi = phase.atom[ia];
            //_p(ia); _parr(atomi.pos,3);
            //_parr(atomi.u,6);

            is = atomi.iscat;

            //_pp(allowed_i[ip][is]); _pp(allowed_j[ip][is]);
            if ( !allowed_i[ip][is] && !allowed_j[ip][is] ) continue;

	    for (sph.rewind(); !sph.finished(); sph.next())
	    {
		for (iatom = ia; iatom != phase.natoms; ++iatom)
		{
		    Atom& atomj = phase.atom[iatom];

		    for (jj=0; jj<3; jj++)
		    {
			dd[jj] = atomi.pos[jj] - atomj.pos[jj] -
				 sph.mno[jj]*phase.icc[jj];
			d[jj] = dd[jj] * phase.a0[jj];
		    }
		    dist2 = phase.skalpro(dd,dd);
		    //cout << dd[0] << " " << dd[1] << " " << dd[2] << " " << dist2 << " " << rmin2 << " " << rmax2 <<endl;
		    //cout << dist2 << " " << rmin2 << " " << rmax2 <<endl;
		    if ( (dist2 >= rmin2) && (dist2 <= rmax2) )
		    {
			//if ((!k) && (j==-1) && (i==1) ) {_pp(dist2); _pp(ja); }

			js = atomj.iscat;

			//_pp(allowed_i[is]); _pp(allowed_j[js]); _pp(allowed_j[is]); _pp(allowed_i[js]);
			if ( (allowed_i[ip][is] &&  allowed_j[ip][js] )
				|| (allowed_j[ip][is] && allowed_i[ip][js]) )
			{

			    //------ Setting up 'thermal' Gaussian
			    sigmap = sqrt(phase.msdAtoms(atomi, atomj, dd));
			    // neglect unphysical summed square displacements
			    if (sigmap <= 0.0)	continue;

			    /* // old crappy code 
			    sigma02 =
				(atomi.u[0]+atomj.u[0])*sqr(d[0]) +
				(atomi.u[1]+atomj.u[1])*sqr(d[1]) +
				(atomi.u[2]+atomj.u[2])*sqr(d[2]) +
				(atomi.u[3]+atomj.u[3])*d[0]*d[1]*2.0 +
				(atomi.u[4]+atomj.u[4])*d[0]*d[2]*2.0 +
				(atomi.u[5]+atomj.u[5])*d[1]*d[2]*2.0;

			    if (sigma02 <= 0)
				continue;   // neglect contribution
			    else
				sigmap = sqrt(sigma02/dist2);
			    */

			    // dist is distance r_ij
			    dist = sqrt(dist2);

			    //- PDF peak width modifications - new 07/2003

#if defined(NEW_SHARP)
			    // Computation of peak sharpening sigma
			    // sigma = sigmap* sqrt(1 - delta/r_ij^2 - gamma/r_ij + qalp*r_ij^2)
			    double corfact;
			    corfact = 1 - phase.delta/dist2
				- phase.gamma/dist + sqr(qalp)*dist2;

			    // if sigma negative: set it back to 1e-5 just for this point
			    // note: derivative will be incorrect in this case
			    if (corfact <= 0)
			    {
				continue;    // neglect contribution
			    }
			    else
			    {
				sigma = sigmap * sqrt(corfact);
			    }
#else
			    // Computation of peak sharpening sigma
			    // sigma = sqrt(sqr(sigmap) - delta/r_ij^2 - gamma/r_ij + qalp*r_ij^2)
			    double sigma2;
			    sigma2 = sqr(sigmap) - phase.delta/dist2
				- phase.gamma/dist + sqr(qalp)*dist2;

			    // if sigma negative: set it back to 1e-5 just for this point
			    // note: derivative will be incorrect in this case
			    if (sigma2 <= 1e-10)
			    {
				sigma = 1e-5;    // neglect contribution
			    }
			    else
			    {
				sigma = sqrt(sigma2);
			    }
#endif
			    // apply rcut if requested
			    if (dist < phase.rcut)
			    {
				sigma *= phase.srat;
			    }

			    // The gaus curve is computed up to distance of 5 sigma
			    igaus   = 1 + nint(5.0*sigma/deltar);

			    gnorm   = 1.0/(sqrt(zpi)*sigma);
			    ampl    = atomi.occ*atomj.occ*phase.weight[is]*phase.weight[js];
			    if (iatom != ia)	ampl += ampl;

			    rb = max(rcmin,dist-5.0*sigma);
			    re = min(rcmax,dist+5.0*sigma);

			    ib = nint((rb-rmin)/deltar);
			    ie = nint((re-rmin)/deltar);

			    for(ig=ib; ig<=ie; ig++)
			    {
				totcalc++;
				rk = rmin + ig*deltar;
				rg = rk-dist;
				gaus = gnorm * exp(-0.5*sqr(rg/sigma));
				ppp[ig] += ampl*gaus;

				// if derivative are needed
				if (ldiff)
				{
				    pdf_derivative(phase,ia,iatom,rk,
					    sigma,sigmap,dist,d,ampl,gaus,fit,fit_a[ig]);
				}
			    }
			}
		    }
		}
	    }
        }

        //------ - Convert to proper G(r) and add to total PDF

        for (i=ncmin;i<=ncmax;i++)
        {
            r = i*deltar + rmin;

            calc[i][ip] = ppp[i]/phase.np/r - fpi*r*phase.rho0*phase.dnorm;

            if (sigmaq > 0.0)
	    {
                calc[i][ip] *= exp(-sqr(r*sigmaq)/2.0);
	    }

            //if ( (r <= phase.corr_max) || (phase.corr_max <= 0.0) )
            pdftot[i] += skal*phase.skal*calc[i][ip];
        }
    }

    // finalize computation of derivatives if required
    // this HAS to happen before the multiplication by deltar/pi and before
    // the convolution to avoid "double operations" on the derivatives,
    // as this will be taken care of in the derivatives.
    if (ldiff)
    {
        fit_setup_derivatives(fit);  // computes matrix dpdf/dvar
    }

    //for (i=ncmin;i<=ncmax;i++)
    //  cout << i << " " << pdftot[i] << endl;

    // From here on we can restrict ourselves to the range [nfmin,nfmax]
    // for the outerloop

    // Apply Qmax cutoff
    if (qmax > 0.0)
    {
	applyQmaxCutoff(&pdftot[ncmin], ncmax-ncmin+1);
    }

    //------ Subtract reference PDF if required
    //
    if (lref)
    {
        for(i=ncmin; i<=ncmax;i++)   // check limits: f2c
            pdftot[i] -= pdfref[i];
    }

    //_pp(totcalc);
    //_p("Exiting <determine>");
}

/*****************************************************************
    This routine calculates the sums over 'ij' needed for the
    derivatives - this is faster than doing the complete loop
    in 'pdf_determine' again ..
    gaus: gaus[igaus+kk] from <determine>
******************************************************************/
void DataSet::pdf_derivative (Phase &phase, int ia, int ja, double rk, double sigma,
    double sigmap, double dist, double d[3], double ampl,double gaus, Fit &fit,
    double* fit_a_i)
{
    double rd,dg;
    rd = dg = 0.0;
    double s11,s22,s33,s12,s13,s23;
    double drdx, drda, phi, dsdphi, dspdx, dspda;
    int ipar, i, ioffset, joffset;

    //------ Some common calculations

    Atom &atomi=phase.atom[ia], &atomj=phase.atom[ja];

    s11 = atomi.u[0] + atomj.u[0];
    s22 = atomi.u[1] + atomj.u[1];
    s33 = atomi.u[2] + atomj.u[2];
    s12 = atomi.u[3] + atomj.u[3];
    s13 = atomi.u[4] + atomj.u[4];
    s23 = atomi.u[5] + atomj.u[5];

    rd = (rk-dist)/sigma;

    if (dist < phase.rcut)
    {
        phi  = phase.srat;
        dsdphi = sigma/phi;
    }
    else
    {
        phi  = 1.0;
        dsdphi = 0.0;
    }

    // derivative of T_ij wrt sigma
    double T = ampl*gaus;
    double dTds = T/sigma*(rd*rd-1);
    double dTdr = (rk-dist)/sqr(sigma)*T;

#if defined(NEW_SHARP)
    //   define s2 = sp**2 - delta/sqr(rij) - gamma/rij + alpha*sqr(rij)
    //   then s = phi*sqrt(s2)
    double dsds2 = sqr(phi*sigmap)/(2.0*sigma);
    double dsdsp = sigma/sigmap;
    double dsdr = dsds2*(2.0*phase.delta/cube(dist)
                    + phase.gamma/sqr(dist) + 2.0*sqr(qalp)*dist);
#else
    //   define s2 = sp**2 - delta/sqr(rij) - gamma/rij + alpha*sqr(rij)
    //   then s = phi*sqrt(s2)
    double dsds2 = sqr(phi)/(2.0*sigma);
    if (sigma==1e-5) dsds2 = 0;

    double dsdsp = 2.0*sigmap*dsds2;
    double dsdr = dsds2*(2.0*phase.delta/cube(dist)
                    + phase.gamma/sqr(dist) + 2.0*sqr(qalp)*dist);
#endif
    double dspdr = - sigmap/dist;



    /*----------------------------------------------------------------
    //------ Derivatives per atom : x,y,z,u,o
    //- -------------------------------------------------------------

      dTdx = dTds.(dsdsp.(dspdr.drdx+dspdx) + dsdr.drdx) + dTdr.drdx
     */

    ioffset = atomi.offset;
    joffset = atomj.offset;

    // ----- d/dx[ip,ia]

    if ( (fit.refvar[ioffset] != -1) || (fit.refvar[joffset] != -1) )
    {
        drdx = phase.a0[0]/dist*(d[0] + phase.cosg*d[1] + phase.cosb*d[2]);
        dspdx = phase.a0[0]*(d[0]*s11 + d[1]*s12 + d[2]*s13)/sqr(dist)/sigmap;

        dg = dTds*(dsdsp*(dspdr*drdx+dspdx) + dsdr*drdx) + dTdr*drdx;
    }

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] -= dg;
    }


    // ----- d/dy[ip,ia]

    if ( (fit.refvar[ioffset] != -1) || (fit.refvar[joffset] != -1) )
    {
        drdx = phase.a0[1]/dist*(phase.cosg*d[0] + d[1] + phase.cosa*d[2]);
        dspdx = phase.a0[1]*(d[0]*s12 + d[1]*s22 + d[2]*s23)/sqr(dist)/sigmap;

        dg = dTds*(dsdsp*(dspdr*drdx+dspdx) + dsdr*drdx) + dTdr*drdx;
    }

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] -= dg;
    }


    // ----- d/dz[ip,ia]

    if ( (fit.refvar[ioffset] != -1) || (fit.refvar[joffset] != -1) )
    {
        drdx = phase.a0[2]/dist*(phase.cosb*d[0] + phase.cosa*d[1] + d[2]);
        dspdx = phase.a0[2]*(d[0]*s13 + d[1]*s23 + d[2]*s33)/sqr(dist)/sigmap;

        dg = dTds*(dsdsp*(dspdr*drdx+dspdx) + dsdr*drdx) + dTdr*drdx;
    }

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] -= dg;
    }


    // ----- d/du[j,ip,ia] (j=11,22,33,12,13,23)

    /*
        dTdu = dTds.dsdsp.dspdu
    */

    // u[0], u[1], u[2]

    double dspdu = 1/(2.0*sigmap*sqr(dist));

    double fsimp = dTds*dsdsp*dspdu;

    for (i=0; i<3; i++)
    {
        dg = fsimp*sqr(d[i]);

        if ( (ipar=fit.refvar[ioffset++]) != -1)
        {
            fit_a_i[ipar] += dg;
            //cout << i << " " << dg << endl; _parr(d,3);
        }
        //cout << phase.a0[i] << " " << d[i] << " " << ig << " " << i << " " << fit_a_i[ipar] << endl;

        if ( (ipar=fit.refvar[joffset++]) != -1)
	{
            fit_a_i[ipar] += dg;
	}
    }

    // u[3]

    dg = 2.0*fsimp*d[0]*d[1];

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
        //cout << 3 << " " << dg << endl; _parr(d,3);
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    // u[4]

    dg = 2.0*fsimp*d[0]*d[2];

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] += dg;
    }

    // u[5]

    dg = 2.0*fsimp*d[1]*d[2];

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
      fit_a_i[ipar] += dg;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
      fit_a_i[ipar] += dg;
    }


    // ----- d/d occ[ip,ia], d/d occ[ip,ja]

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        fit_a_i[ipar] += T/atomi.occ;
    }

    if ( (ipar=fit.refvar[joffset++]) != -1)
    {
        fit_a_i[ipar] += T/atomj.occ;
    }

    //------ ----------------------------------------------------------------
    //------ Derivatives per phase : lat,delta,gamma
    //------ ----------------------------------------------------------------

    // ----- d/d(lattice parameter a)

    //_pp(&phase.a0[0]); _pp(fit.var[0]); _pp(fit.fit(phase.a0[0]));
    ioffset = phase.offset;

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = (d[0] + phase.cosg*d[1] + phase.cosb*d[2])*d[0]/phase.a0[0]/dist;
        dspda = (d[0]*s11 + d[1]*s12 + d[2]*s13)*d[0]/phase.a0[0]/sqr(dist)/sigmap;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda + dTds*dsdsp*dspda;

        fit_a_i[ipar] += dg;
    }

    // ----- d/d(lattice parameter b)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = (phase.cosg*d[0] + d[1] + phase.cosa*d[2])*d[1]/phase.a0[1]/dist;
        dspda = (d[0]*s12 + d[1]*s22 + d[2]*s23)*d[1]/phase.a0[1]/sqr(dist)/sigmap;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda + dTds*dsdsp*dspda;

        fit_a_i[ipar] += dg;
    }

    // ----- d/d(lattice parameter c)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = (phase.cosb*d[0] + phase.cosa*d[1] + d[2])*d[2]/phase.a0[2]/dist;
        dspda = (d[0]*s13 + d[1]*s23 + d[2]*s33)*d[2]/phase.a0[2]/sqr(dist)/sigmap;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda + dTds*dsdsp*dspda;

        fit_a_i[ipar] += dg;
    }


    // ----- d/d(lattice angle alpha)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = - rad*phase.sina*d[1]*d[2]/dist;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda;

        fit_a_i[ipar] += dg;
    }

    // ----- d/d(lattice angle beta)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = - rad*phase.sinb*d[0]*d[2]/dist;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda;

        fit_a_i[ipar] += dg;
    }

    // ----- d/d(lattice angle gamma)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        drda = - rad*phase.sing*d[0]*d[1]/dist;

        dg = (dTds*(dsdsp*dspdr+dsdr) + dTdr)*drda;

        fit_a_i[ipar] += dg;
    }

    // Derivatives wrt the peak sharpening parameters

    // ----- d/d(delta)


    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = -dTds*dsds2/sqr(dist);
        fit_a_i[ipar] += dg;
    }

    // ----- d/d(gamma)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = -dTds*dsds2/dist;
        fit_a_i[ipar] += dg;
    }

    // ----- d/d(csca)

    ioffset++;


    // ----- d/d(srat)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = dTds*dsdphi;
        fit_a_i[ipar] += dg;
    }


    //------ ----------------------------------------------------------------
    //------ Derivatives per data set : dsca, qsig, qalp
    //------ ----------------------------------------------------------------

    ioffset = this->offset;

    // ----- d/d(dsca[is])

    ioffset++;

    // ----- d/d(qsig[is])

    ioffset++;

    // ----- d/d(qalp[ip])

    if( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = dTds*dsds2*(2.0*qalp*sqr(dist));
        fit_a_i[ipar] += dg;
    }

}


/****************************************************************
    Setup for various arrays and functions for PDF calculation.
       for phase np and data set ns
*****************************************************************/
void Phase::setup_weights(bool lout, bool lxray)
{
    int ia, is;

    matrix<double> scat(nscat,9);
    weight.resize(nscat);


    if (lout) cout << " Setting up PDF segment ...";

    // ------ Setting up weighting (b(i)b(j)/<b**2>)

    dlink(scat, lxray);

    // weights (unnormalized)
    for (is=0; is<nscat; is++)
        weight[is] = scatteringFactor(scat[is], lxray);

    // the function gets computed way too many times!!!  Once for each atom and twice for
    // each atom type
    bave = 0.0;
    for (ia=0; ia<natoms; ia++)
        bave += atom[ia].occ * weight[atom[ia].iscat];

    bave /= np;
    //_pp(bave);

    // normalize the weights
    for (is=0; is<nscat; is++)
    {
        weight[is] /= bave;
        //_pp(weight[is][js]);
    }
}

/************************************************************************
* Extend r-range by 6 ripples of sinc before applying Qmax cutoff. 
* Contributions from peaks outside should be less than 1.5%.
*
* todo: should also extend by the 5*max(sigma)
************************************************************************/
void DataSet::extendCalculationRange(bool lout)
{
    const int nripples = 6;
    // initialize calculation range to fitting range
    rcmin = rfmin;
    rcmax = rfmax;
    ncmin = nint((rcmin - rmin)/deltar);
    ncmax = nint((rcmax - rmin)/deltar);
    // check if Qmax cutoff is applied
    if (!(0.0 < qmax))	return;
    // get extension width
    double rext = nripples*2*M_PI/qmax;
    rcmin = max(rmin, rfmin - rext);
    rcmax = rfmax + rext;
    ncmin = nint((rcmin - rmin)/deltar);
    ncmax = nint((rcmax - rmin)/deltar);
    if (lout)
    {
	cout << " Extending PDF search distance to " <<
	    rcmin << " -> " << rcmax << " A ...\n";
    }
}

/************************************************************************
 * diameter of sphere that can enclose primitive cell
 ************************************************************************/
double Phase::circum_diameter()
{
    const double epsilond = sqrt(numeric_limits<double>().epsilon());
    if (atom.empty())	return 0.0;
    double center[3] = {0.0, 0.0, 0.0};
    for (vector<Atom>::iterator ai = atom.begin(); ai != atom.end(); ++ai)
    {
	for (int i = 0; i != 3; ++i)	center[i] += ai->pos[i];
    }
    for (int i = 0; i != 3; ++i)    center[i] /= natoms;
    double maxd = 0.0;
    for (vector<Atom>::iterator ai = atom.begin(); ai != atom.end(); ++ai)
    {
	double dd[3];
	for (int i = 0; i !=3 ; ++i)	dd[i] = ai->pos[i] - center[i];
	double d = sqrt(skalpro(dd,dd));
	if (d > maxd)	maxd = d;
    }
    maxd = maxd*(1.0+epsilond) + epsilond;
    return 2*maxd;
}

/************************************************************************
 * combined mean square displacements of 2 atoms along lattice vector
 ************************************************************************/
double Phase::msdAtoms(const Atom& ai, const Atom& aj, double* vl)
{
    // msd = transpose(gten*vln) * ar[i]*U[i,j]*ar[j] * gten*vln
    // normalize vl
    double vlnorm = sqrt(skalpro(vl,vl));
    double vln[3] = { vl[0]/vlnorm, vl[1]/vlnorm, vl[2]/vlnorm };
    // combine squared atom displacements
    double U[6];
    for (size_t i = 0; i != 6; ++i)
    {
	U[i] = ai.u[i] + aj.u[i];
    }
    // Un = ar[i]*U[i,j]*ar[j]
    double Un[6] = { U[0]*ar[0]*ar[0], U[1]*ar[1]*ar[1], U[2]*ar[2]*ar[2],
		     U[3]*ar[0]*ar[1], U[4]*ar[0]*ar[2], U[5]*ar[1]*ar[2] };
    // rhs = gten*vln
    double rhs[3] = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i != 3; ++i)
    {
	for (size_t j = 0; j != 3; ++j)
	{
	    rhs[i] += gten[i][j] * vln[j];
	}
    }
    // msd = transpose(rhs) * Un * rhs
    double msd;
    msd = Un[0]*rhs[0]*rhs[0] + Un[1]*rhs[1]*rhs[1] + Un[2]*rhs[2]*rhs[2] +
	2*Un[3]*rhs[0]*rhs[1] + 2*Un[4]*rhs[0]*rhs[2] + 2*Un[5]*rhs[1]*rhs[2];
    return msd;
}

/*********************************
------ - Save fit results
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the
    saved file. Actually saving the file is
    optional.
*********************************/
string PdfFit::save_res(string fname)
{
    //check to see if a refinement has even been done.
    //after a refinement is finished, fit.iter = -1
    if(fit.iter == 0)
    {
        throw unassignedError("Refinement must be performed first");
        return "";
    }
    ofstream fout;
    stringstream outfilestream;
    string outfilestring = "";

    outfilestream << " " << string(78,'=') << endl
    << " PDF REFINEMENT \n"
    << "   Using PDFFIT version : " << version << endl
    << " " << string(78,'=') << endl;

    for(int ip=0; ip<nphase; ip++)
        (*phase[ip]).output(outfilestream);

    for (int is=0; is<nset; is++)
        (*set[is]).output(outfilestream);

    fit.output(outfilestream);

    outfilestream << " " << string(78,'=') << endl;

    if( fname != "" )
    {
        fout.open(fname.c_str());

        if (!fout)
	{
            //warning("save_res: cannot create output file");
            throw IOError("Cannot create output file");
            return outfilestring;
        }

        cout << " Saving fit results to file : " << fname << endl;

        fout << outfilestream.str();
        //fout.setf(ios::showpoint);

        fout.close();
    }
    else
    {
//        cout << " Not saving fit results to file." << endl;
    }

    return outfilestream.str();
}

/*************************************************
 -------  Save PDF, structure or complete result
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the
    saved file. Actually saving the file is
    optional.
**************************************************/
string PdfFit::save_pdf(int iset, string fname)
{
    string outfilestring = "";

    //------ - Save PDF (G(r))

    if ( (iset < 1) || (iset > nset) )
    {
        //warning("save_pdf: data set does not exist\n");
        throw unassignedError("data set does not exist");
        return outfilestring;
    }
    else if( fname != "" )
    {
        ofstream fout;
        fout.open(fname.c_str());

        if (!fout)
	{
            //warning("save_pdf: cannot create output file");
            throw IOError("cannot create output file");
            return outfilestring;
        }

        cout << " Saving PDF data set " << iset << " to file : " << fname << endl;

        outfilestring = set[iset-1]->build_pdf_file();

        fout << outfilestring;
        fout.close();
    }
    else
    {
//        cout << " Not saving PDF data set " << iset << " to file." << endl;
        outfilestring = set[iset-1]->build_pdf_file();
    }

    return outfilestring;
}


/*
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the
    saved file.
*/
string DataSet::build_pdf_file()
{
    string blank=string(4,' ');
    stringstream outfilestream;
    outfilestream.setf(ios::showpoint);

    for (int i=nfmin;i<=nfmax;i++)
    {
        double r = i*deltar + rmin;
        outfilestream << setw(12) << r << blank << setw(12) << pdftot[i] << blank << setw(12) << 0.0
                      << blank << setw(12) << 1.0/sqrt(wic[i]) << blank << setw(12)
                      << obs[i]-pdftot[i] << endl;
    }

    return outfilestream.str();
}

//------ - Save DIF file (Gobs-G(r))
/*
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the
    saved file. Actually saving the file
    is optional.
*/

string PdfFit::save_dif(int iset, string fname)
{

    ofstream fout;
    string outfilestring = "";

    if ( (iset < 1) || (iset > nset) )
    {
        //warning("save_dif: data set does not exist\n");
        throw unassignedError("Data set does not exist");
    }
    else if (fname != "" )
    {
        fout.open(fname.c_str());

        if (!fout)
	{
            //warning("save_dif: cannot create output file");
            throw IOError("Cannot create output file");
            return outfilestring;
        }

        cout << " Saving difference data set " << iset << " to file : " << fname << endl;

        outfilestring = set[iset-1]->build_dif_file();

        fout << outfilestring;
        fout.close();
    }
    else
    {
//        cout << " Not saving difference data set " << iset << " to file." << endl;
        outfilestring = set[iset-1]->build_dif_file();
    }

    return outfilestring;
}

/*
    Thu Oct 13 2005 - CLF
    Changed code to return a string of the saved file.
*/
string DataSet::build_dif_file()
{
    string blank=string(4,' ');
    stringstream outfilestream;
    outfilestream.setf(ios::showpoint);

    for (int i=nfmin;i<=nfmax;i++)
    {
        double r = i*deltar + rmin;
        outfilestream << setw(12) << r << blank << setw(12) << obs[i]-pdftot[i] << endl;
    }

    return outfilestream.str();
}

void PdfFit::selphase(int ip)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }

    if (ip==ALL)
    {
        for (int i=0; i<nphase; i++)
            curset->selphase(i, phase[i]);
    }
    else
    {
        if (ip <= nphase)
	{
            curset->selphase(ip-1, phase[ip-1]);
	}
        else
	{
            stringstream eout;
            eout << "Phase " << ip << " undefined";
            throw unassignedError(eout.str());
        }
    }
}

void DataSet::selphase(int ip, Phase *phase)
{
    if ( (int) psel.size() <= ip)
    {
        psel.resize(ip+1);
        allowed_i.resize(ip+1);
        allowed_j.resize(ip+1);
    }

    psel[ip] = phase;

    allowed_i[ip].resize(phase->nscat);
    for (int i=0; i<phase->nscat; i++)
        allowed_i[ip][i] = true;

    allowed_j[ip].resize(phase->nscat);
    for (int i=0; i<phase->nscat; i++)
        allowed_j[ip][i] = true;
}

void PdfFit::pdesel(int ip)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }

    if (ip==ALL)
    {
        for (int i=0; i<nphase; i++)
            curset->psel[ip] = NULL;
    }
    else
    {
        if (ip <= nphase)
	{
            curset->psel[ip-1] = NULL;
	}
        else
	{
            stringstream eout;
            eout << "phase " << ip << " undefined";
            throw unassignedError(eout.str());
        }
    }
}


/********************************************
  This routine executes the select command
*********************************************/
void DataSet::selatom(int ip, int i, vector<vector<bool> > &allowed, bool choice)
{
    stringstream eout; //for errors
    if (ip==ALL)
    {
        for(unsigned int j=0; j < psel.size(); j++)
        {
            if(i==ALL)
            {
                for(int k=0; k<psel[j]->nscat; k++)
                    allowed[j][k] = choice;
            }
            else
            {
                //cout << "Warning: phase = ALL requires atom = ALL" << endl;
                throw ValueError("phase = ALL requires atom = ALL");
                return;
            }
        }
    }
    else if ( (ip>=0) && (ip< (int) psel.size()) )
    {
        if(i==ALL)
        {
            for(int k=0; k<psel[ip]->nscat; k++)
                allowed[ip][k] = choice;
        }
        else if ( (i>=0) && (i<psel[ip]->nscat) )
	{
            allowed[ip][i] = choice;
	}
        else
	{
            eout << "atom type " << i+1 << " undefined\n";
            throw unassignedError(eout.str());
	}
    }
    else
    {
        eout << "phase " << ip+1 << " undefined or not selected\n";
        throw unassignedError(eout.str());
    }
}

void PdfFit::isel(int ip, int i)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }
    else
    {
        curset->selatom((ip==ALL)?ALL:ip-1, (i==ALL)?ALL:i-1, curset->allowed_i, true);
    }
}

void PdfFit::idesel(int ip, int i)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }
    else
    {
        curset->selatom((ip==ALL)?ALL:ip-1, (i==ALL)?ALL:i-1, curset->allowed_i, false);
    }
}

void PdfFit::jsel(int ip, int i)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }
    else
    {
        curset->selatom((ip==ALL)?ALL:ip-1, (i==ALL)?ALL:i-1, curset->allowed_j, true);
    }
}

void PdfFit::jdesel(int ip, int i)
{
    if (!curset)
    {
        //warning("No data set selected");
        throw unassignedError("No data set selected");
        return;
    }
    else
    {
        curset->selatom((ip==ALL)?ALL:ip-1, (i==ALL)?ALL:i-1, curset->allowed_j, false);
    }
}


/*****************************************
    Wed Oct 12 2005 - CLF
    Reads observed PDF from arrays.
******************************************/
int PdfFit::read_data_arrays(Sctp t, double qmax, double sigmaq,
        int length, double * r_data, double * Gr_data, double * dGr_data, string _name)
{
    DataSet *_set = new DataSet, &set=*_set;

    set.read_data_arrays(nset+1, t, qmax, sigmaq, false, length, r_data, Gr_data, dGr_data, _name);

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        set.selphase(ip, this->phase[ip]);

    this->set.push_back(&set);
    nset++;
    setdata(nset);

    return 1;
}

/*****************************************
    Wed Oct 12 2005 - CLF
    Reads observed PDF from a c-style string.
******************************************/

int PdfFit::read_data_string(string& buffer, Sctp t, double qmax, double sigmaq, string _name)
{
    DataSet *_set = new DataSet, &set=*_set;

    try{
        set.read_data_string(nset+1, buffer, t, qmax, sigmaq, false);
    }
    catch(Exception e) {
        delete &set;
    //  cout << "Error reading data <" << _name << ">: "
    //      << e.GetMsg() << " --> no new data set allocated\n";
        throw;
        return 0;
    }

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        set.selphase(ip, this->phase[ip]);

    this->set.push_back(&set);
    nset++;
    setdata(nset);

    return 1;
}

/*****************************************
    Reads observed PDF as xy ASCII file.
******************************************/

int PdfFit::read_data(string datafile, Sctp t, double qmax, double sigmaq)
{
    DataSet *_set = new DataSet, &set=*_set;

    try{
        set.read_data(nset+1, datafile, t, qmax, sigmaq, false);
    }
    catch(Exception e) {
        delete &set;
    //  cout << "Error in data file <" << datafile << ">: "
    //      << e.GetMsg() << " --> no new data set allocated\n";
        throw;
        return 0;
    }

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        set.selphase(ip, this->phase[ip]);

    this->set.push_back(&set);
    nset++;
    setdata(nset);

    return 1;
}

static void extract_key(string line, string key);
vector<double> res_para;

/* Wed Oct 12 2005 - CLF
 * Using read_data_arrays adds functionality
 * to pdffit2, allowing one to read data that is alread stored as arrays.
 */
void DataSet::read_data_arrays(int _iset, Sctp t, double _qmax, double _sigmaq,
        bool lref, int length, double * r_data, double * Gr_data,
        double * dGr_data, string _name )
{
    iset = _iset;

    //------ Now analyse given parameters

    //------ - get exp. method (neutron/x-ray)

    if (t == N)
    {
        lxray = false;
    }
    else if (t == X)
    {
        lxray = true;
    }

    //------ - get QMAX and QSIG

    qmax = _qmax;
    sigmaq = _sigmaq;

    //------ - read possible reference PDF filename

    if (lref)
    {
        this->lref = true;
        //call del_params(1,ianz,cpara,lpara,maxw)
        //call do_build_name(ianz,cpara,lpara,werte,maxw,1)
        //rfile  = cpara(1)
        //lrfile = lpara(1)
    }
    else
    {
        this->lref = false;
    }

    //------ Finally we actually read the data

    bool lwei = true;

    /* Only really care about G(r) and dG(r), the uncertainty */
    for( int i = 0; i < length; i++ )
    {
        double wic;

        if( dGr_data == NULL )
        {
            wic = 1.0;
            lwei = false;
        }
        else
        {
            wic = 1.0/sqrt( dGr_data[i] );
        }

        this->obs.push_back( Gr_data[i] );
        this->wic.push_back(wic);

    }

    cout << " Reading data from arrays...\n";

    //  FIRST COLUMN IS SEEMINGLY NOT USED
    // WHAT HAPPENS IF THE DATA FILE IS NOT EQUIDISTANT???
    // WHAT HAPPENS IF THE DATA FILE IS NOT SORTED CORRECTLY???
    rmin = rfmin = r_data[0];
    rmax = rfmax = r_data[length - 1];
    bin = length;
    deltar = (rmin-rmax)/double(bin-1);
    name = _name;

    cout << " Read PDF data set " << iset << "  (r = " << rmin << " to " << rmax << " A, " << bin << " points) ...\n";
    if (!lwei)
    {
        cout << " No sigmas for G(r) found, using unit weights ...\n";
    }
    cout << endl;

    return;
}

/* Wed Oct 12 2005 - CLF
 * Using read_data_string adds functionality
 * to pdffit2, allowing one to read data that has already been loaded.
 *
 *
 *  Thu Nov  3 2005 - CLF
 *  Need to add some data checking routines!
 */
void DataSet::read_data_string(int _iset, string& buffer, Sctp t, double _qmax,
        double _sigmaq, bool lref, string _name)
{
    //ifstream fdata;
    istringstream fdata( buffer );

    string line;
    int ncol = 2;

    iset = _iset;

    //------ Now analyse given parameters

    //------ - get exp. method (neutron/x-ray)

    if (t == N)
    {
        lxray = false;
    }
    else if (t == X)
    {
        lxray = true;
    }

    //------ - get QMAX and QSIG

    qmax = _qmax;
    sigmaq = _sigmaq;

    //------ - read possible reference PDF filename

    if (lref)
    {
        this->lref = true;
        //call del_params(1,ianz,cpara,lpara,maxw)
        //call do_build_name(ianz,cpara,lpara,werte,maxw,1)
        //rfile  = cpara(1)
        //lrfile = lpara(1)
    }
    else
    {
        this->lref = false;
    }

    //------ Read observed PDF for given plane

    //fdata.open(pfile.c_str());

    //if (!fdata) throw Exception("file does not exist");

    //if (lref) call oeffne (18,rfile,'old',.false.)

    //------ Has the file a prepended History ?
    //------ Get information from history and store in res[i]

    //------ res[1] = Temperature

    // reset res_para vector
    res_para.clear();

    getline(fdata, line);

    if (!line.compare(0,7,"History"))
    {
        cout << " History information found ...\n";

        while(1)
        {
            getline(fdata, line);  //   if err 999
            extract_key(line,"temp=");
            extract_key(line,"Qmax=");
            if (strcmp(line,"##### start data",16)) break;
        }
    }

    //------ Any other header lines starting with # ?

    while (line[0] == '#')
    {
        getline(fdata, line);
    }

    //------ Finally we actually read the data

    double re = 0.0;
    double ra = 0.0;
    bool lwei = true;

    while (1)
    {
        double val, obs, wic;
        istringstream sline(line);

        sline >> re >> obs;

        ncol = 2;

        if (!this->obs.size()) ra = re;

        sline >> val;

        // try to read a 3rd column
        if (sline)
        {
            ncol++;

            if  (val > 0.0)
	    {
                wic = 1.0/sqr(val);
	    }

            // try to read a 4th column in the line
            sline >> val;

            // if 4 columns are present: reshuffle variables (dy) -> (dx,dy)
            if (sline && (val > 0.0) )
            {
                ncol++;
                wic = 1.0/sqr(val);
            }
        }
        else
        {
            wic = 1.0;
            lwei = false;
        }

        //if (lref) then
        //read (18,*,end=20,err=999) de,pdf_ref(ip,iset)
        //if (abs(re-de).gt.1e-5) goto 9999
        //endif

        this->obs.push_back(obs);
        this->wic.push_back(wic);

        if (not getline(fdata, line))	break;
    }

    cout << " Reading " << ncol << " columns ...\n";

    //if (lref) fref.close();


    //  FIRST COLUMN IS SEEMINGLY NOT USED
    // WHAT HAPPENS IF THE DATA FILE IS NOT EQUIDISTANT???
    // WHAT HAPPENS IF THE DATA FILE IS NOT SORTED CORRECTLY???
    rmin = rfmin = ra;
    rmax = rfmax = re;
    bin = this->obs.size();
    deltar = (re-ra)/double(bin-1);
    name = _name;

    //if (lref) pdf_rname(iset)  = rfile(1:lrfile)

    //if (lref) then
    //  cout << " Read PDF difference data set " << iset << "  (r = "
    //  << ra << " to " << re << " A, " << bin << " points) ...\n";
    //else
      cout << " Read PDF data set " << iset << "  (r = " << ra
            << " to " << re << " A, " << bin << " points) ...\n";
    //endif

    if (!lwei) cout << " No sigmas for G(r) found, using unit weights ...\n";

    cout << endl;

    //for (int i=0; i<bin; i++)
    //  _pp(obs[i]);

    return;
}

/* Wed Oct 12 2005 - CLF
 * Using read_data and the above read_data_string adds functionality
 * to pdffit2, allowing one to read data that has already been loaded.
 */
void DataSet::read_data(int _iset, string pfile, Sctp t, double _qmax,
        double _sigmaq, bool lref)
{
    ifstream fdata;

    //------ Read observed PDF for given plane
    fdata.open(pfile.c_str());
    if (!fdata) throw IOError("File does not exist");

    // Read the file into a buffer and send it off to read_data_string function
    //char * cbuffer;
    int length;
    fdata.seekg(0, ios::end);
    length = fdata.tellg();
    fdata.seekg(0, ios::beg);
    //cbuffer = new char [length];
    char cbuffer[length];
    fdata.read( cbuffer, length );
    fdata.close();
    string buffer(cbuffer, length);

    read_data_string( _iset, buffer, t, _qmax, _sigmaq, lref, pfile );

    return;

}

/************************************************
    Gets numbers from history part of data file
*************************************************/
static void extract_key(string line, string key)
{
    string::size_type is;

    is = line.find(key);

    if (is != string::npos)
    {
        istringstream sline(line.substr(is+key.size(),line.size()));

        res_para.push_back(dget(sline));
    }
}


/*********************************
    Sets R-range for fitting
**********************************/
void PdfFit::range(int is, double rmin, double rmax)
{
    if( rmin >= rmax )
    {
	throw ValueError("rmin must be < rmax");
    }
    if (is == ALL)
    {
	for(is=0;is<nset;is++)
	    (*set[is]).range(rmin,rmax);
    }
    else
    {
	if ( (is>=1) && (is<=nset) )
	{
	    (*set[is-1]).range(rmin,rmax);
	}
	else
	{
	    //cout << "Invalid data set number\n";
	    throw ValueError("Invalid data set number");
	}
    }
}

void DataSet::range(double rmin, double rmax)
{
    DataSet &set=*this;

    if ( (rmin >= set.rmin) && (rmin <= set.rmax)
        && (rmax <= set.rmax) && (rmin < rmax) )
    {
        set.rfmin = rmin;
        set.rfmax = rmax;
    }
    else
    {
        throw ValueError("Range outside data set limits");
    }
}


#if defined(FORTRAN)
c*****7*****************************************************************
	subroutine do_xray (zeile,lp)
c-
c	Setting Q to calculate f(xray) ..
c+
	implicit      	none
c
	include 	'config.inc'
	include 	'prompt.inc'
	include 	'pdf.inc'
	include 	'errlist.inc'
c
	integer       	maxw
	parameter    	(maxw=2)
c
	character*(*)	zeile
	integer		lp
c
	character*200 	cpara(maxw)
	integer       	lpara(maxw),length
	integer       	ianz
	real          	werte(maxw)
c
	call get_params (zeile,ianz,cpara,lpara,maxw,lp)
	if (ier_num.ne.0) return
c
	if     (ianz.eq.0) then
	  write (output_io,1000) pdf_xq
	elseif (ianz.eq.1) then
	  call ber_params(ianz,cpara,lpara,werte,maxw)
	  if (ier_num.ne.0) return
	  pdf_xq = werte(1)
	else
	  ier_num = -6
	  ier_typ = ER_COMM
	endif
c
1000	format (1x,'------ > Current Q for X-ray weights : ',g14.6)
	end
#endif

