/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Mixed methods for PDF calculation from PdfFit, DataSet and Phase
*
* Comments: Up to date with 1.3.10 Fortran version.
*	    What a spagetti.
*
***********************************************************************/

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <valarray>
#include <cassert>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "PointsInSphere.h"
#include "StringUtils.h"
#include "LocalPeriodicTable.h"
#include "MathUtils.h"
#include "ShapeFactors.h"

#include "pdffit.h"

using NS_PDFFIT2::pout;

#define NEW_SHARP

/*****************************************************************
  Allocating space for dummy dataset when calculating PDF
  without data
******************************************************************/
void PdfFit::alloc(char tp, double qmax, double qdamp, double rmin, double rmax, int bin)
{

    DataSet* pds = new DataSet(this);
    int i;
    if (rmax < rmin || rmin < 0 || rmax < 0)
    {
        throw ValueError("Check rmin, rmax");
    }
    if (bin <= 1)
    {
        throw ValueError("bin must be > 1");
    }
    if( qmax < 0.0 )
    {
        throw ValueError("qmax must be >= 0");
    }
    if( qdamp < 0 )
    {
        throw ValueError("qdamp must be >= 0");
    }

    // 2009-03-11 PJ: Removed check if a structure has been loaded before.
    // It can be now loaded after calling alloc.

    pds->iset = nset+1;

    pds->scattering_type = tp;

    pds->qmax   = qmax;
    pds->qdamp = qdamp;
    pds->rmin   = pds->rfmin = rmin;
    pds->rmax   = pds->rfmax = rmax;
    pds->bin    = bin;
    pds->deltar = (rmax-rmin)/double(bin-1);
    pds->name   = "Dummy set";

    pds->obs.resize(bin);
    pds->wic.resize(bin);

    for (i=0; i<bin; i++)
    {
        pds->obs[i] = 0.0;
        pds->wic[i] = 1.0;
    }

    *pout << " Allocated PDF data set " << pds->iset << "  (r = "
        << rmin << " to " << rmax << " A, " << bin << " points) ..." << endl;

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        pds->selphase(ip, this->phase[ip]);

    this->datasets.push_back(pds);
    nset++;
    setdata(nset);
}

// cut-away high Q harmonics using fast Fourier transformation
void DataSet::applyQmaxCutoff(double* y, size_t len)
{
    // pad y with the same number of zeros up to the next power of 2
    size_t padlen = 2*len;
    // MS compatibility fix
    padlen = (size_t)  pow(2.0, ceil(log2(padlen))) ;
    // ycpad is complex, so it needs to be twice as long
    valarray<double> ycpad_array(2*padlen);
    double* ycpad = &(ycpad_array[0]);
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
    //int kk, ibin, ip;
    int ib, ie, ig;
    double rmax2, rmin2, gaus, rk, rb,re, rtot, ract;
    vector<double> ppp;
    double dd[3], d[3], dist2, dist, rg, r;
    double sigmap, sigma, gnorm, ampl;
    long totcalc=0;

    if (lout) *pout << " Calculating PDF ...\n";

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

        Phase& phase = *psel[ip];

        phase.setup_weights(scattering_type);

	// build flag arrays for ignored i and j atoms.  The skip expression
	// is symmetric in ignore_i, ignore_j thus it is OK to do summation
	// only over the upper triangular part of the pair matrix.
	valarray<bool> ignore_i(phase.natoms);
	valarray<bool> ignore_j(phase.natoms);
	for (int aidx = 0; aidx < phase.natoms; ++aidx)
	{
	    ignore_i[aidx] = phase_ignore_i[&phase].count(aidx);
	    ignore_j[aidx] = phase_ignore_j[&phase].count(aidx);
	}

        // ------ Get the ratio total pairs/selected weight in structure

        rtot = 0.0;
        ract = 0.0;
	for (int iidx = 0; iidx < phase.natoms; ++iidx)
	{
	    Atom& ai = phase.atom[iidx];
	    for (int jidx = iidx; jidx < phase.natoms; ++jidx)
	    {
		Atom& aj = phase.atom[jidx];
		double halfloopscale = (iidx == jidx) ? 1.0 : 2.0;
                rtot += halfloopscale * ai.weight * aj.weight;
		// skip ignored pair
		bool skip = (ignore_i[iidx] || ignore_j[jidx]) &&
		    (ignore_i[jidx] || ignore_j[iidx]);
		if (skip)   continue;
		ract += halfloopscale * ai.weight * aj.weight;
	    }
	}
        phase.dnorm = ract/rtot;

        //--------------------------------------------------------------

	fill(ppp.begin() + ncmin, ppp.begin() + ncmax + 1, 0.0);

	// calculate range for PointsInSphere sequencer
	// (negative rsphmin is no problem)
	double buffzone = phase.circum_diameter();
	double rsphmin = sqrt(rmin2) - buffzone;
	double rsphmax = sqrt(rmax2) + buffzone;
        // limit rsphmax, when stepcut has been set for the phase
        if (phase.stepcut > 0.0)
        {
            rsphmax = min(rsphmax, phase.stepcut + buffzone);
        }
	PointsInSphere sph( rsphmin, rsphmax, phase.a0[0]*phase.icc[0],
		phase.a0[1]*phase.icc[1], phase.a0[2]*phase.icc[2],
		phase.win[0], phase.win[1], phase.win[2] );

	// loop only over selected atom indexes
	for (int iidx = 0; iidx < phase.natoms; ++iidx)
	{
	    Atom& ai = phase.atom[iidx];
	    for (int jidx = iidx; jidx < phase.natoms; ++jidx)
	    {
		Atom& aj = phase.atom[jidx];

		// skip ignored pair
		bool skip = (ignore_i[iidx] || ignore_j[jidx]) &&
		    (ignore_i[jidx] || ignore_j[iidx]);
		if (skip)   continue;

		for (sph.rewind(); !sph.finished(); sph.next())
		{
		    for (int i=0; i<3; i++)
		    {
			dd[i] = ai.pos[i] - aj.pos[i] -
				 sph.mno[i]*phase.icc[i];
			d[i] = dd[i] * phase.a0[i];
		    }
		    dist2 = phase.skalpro(dd,dd);

                    // check if pair distance is within dataset r limits
                    if ((dist2 < rmin2) || (dist2 > rmax2))    continue;

                    // dist is distance r_ij
                    dist = sqrt(dist2);

                    //------ Setting up 'thermal' Gaussian
                    sigmap = sqrt(phase.msdAtoms(ai, aj, dd));
                    // neglect unphysical summed square displacements
                    if (sigmap <= 0.0)	continue;

                    //- PDF peak width modifications - new 07/2003

#if defined(NEW_SHARP)
                    // Computation of peak sharpening sigma
                    // sigma = sigmap* sqrt(1 - delta2/r_ij^2 - delta1/r_ij + qbroad*r_ij^2)
                    double corfact;
                    corfact = 1 - phase.delta2/dist2
                        - phase.delta1/dist + sqr(qbroad)*dist2;

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
                    // sigma = sqrt(sqr(sigmap) - delta2/r_ij^2 - delta1/r_ij + qbroad*r_ij^2)
                    double sigma2;
                    sigma2 = sqr(sigmap) - phase.delta2/dist2
                        - phase.delta1/dist + sqr(qbroad)*dist2;

                    // if sigma2 negative: set it back to 1e-5 just for this point
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
                        sigma *= phase.sratio;
                    }

                    // The gaus curve is computed up to distance of 5 sigma
                    gnorm   = 1.0/(sqrt(2.0*M_PI)*sigma);
                    ampl    = ai.occ * ai.weight * aj.occ * aj.weight;

                    if (iidx != jidx)   ampl += ampl;

                    rb = max(rcmin, dist - 5.0*sigma);
                    re = min(rcmax, dist + 5.0*sigma);
                    if (phase.stepcut > 0.0)
                    {
                        re = min(re, phase.stepcut);
                    }

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
                            pdf_derivative(phase, ai, aj, rk,
                                    sigma, sigmap, dist, d, ampl, gaus,
                                    fit, fit_a[ig]);
                        }
                    }
		}
	    }
        }

        //------ - Convert to proper G(r) and add to total PDF

        for (int i = ncmin; i<=ncmax; i++)
        {
            r = i*deltar + rmin;

            calc[i][ip] = ppp[i]/phase.np/r - 4.0*M_PI*r*phase.rho0*phase.dnorm;

	    // Q-resolution envelope
            if (qdamp > 0.0)
	    {
                calc[i][ip] *= exp(-sqr(r*qdamp)/2.0);
	    }

	    // PDF envelope for spherical nano particles
            if (phase.spdiameter > 0.0)
            {
                calc[i][ip] *= sphereEnvelope(r, phase.spdiameter);
            }

            // Empirical shape step cutoff of the PDF
            if (phase.stepcut > 0.0)
            {
                double stepfactor = (r <= phase.stepcut) ? 1.0 : 0.0;
                calc[i][ip] *= stepfactor;
            }

            //if ( (r <= phase.corr_max) || (phase.corr_max <= 0.0) )
            pdftot[i] += this->dscale * phase.pscale * calc[i][ip];
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
    //  *pout << i << " " << pdftot[i] << endl;

    // From here on we can restrict ourselves to the range [nfmin,nfmax]
    // for the outerloop

    // Apply Qmax cutoff
    if (qmax > 0.0)
    {
	applyQmaxCutoff(&pdftot[ncmin], ncmax-ncmin+1);
    }
}

/*****************************************************************
    This routine calculates the sums over 'ij' needed for the
    derivatives - this is faster than doing the complete loop
    in 'pdf_determine' again ..
    gaus: gaus[igaus+kk] from <determine>
******************************************************************/
void DataSet::pdf_derivative (Phase &phase,
	const Atom& atomi, const Atom& atomj, double rk, double sigma,
	double sigmap, double dist, double d[3], double ampl,double gaus,
	Fit &fit, double* fit_a_i)
{
    double rd,dg;
    rd = dg = 0.0;
    double s11,s22,s33,s12,s13,s23;
    double drdx, drda, phi, dsdphi, dspdx, dspda;
    int ipar, i, ioffset, joffset;

    //------ Some common calculations

    s11 = atomi.u[0] + atomj.u[0];
    s22 = atomi.u[1] + atomj.u[1];
    s33 = atomi.u[2] + atomj.u[2];
    s12 = atomi.u[3] + atomj.u[3];
    s13 = atomi.u[4] + atomj.u[4];
    s23 = atomi.u[5] + atomj.u[5];

    rd = (rk-dist)/sigma;

    if (dist < phase.rcut)
    {
        phi  = phase.sratio;
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
    //   define s2 = sp**2 - delta2/sqr(rij) - delta1/rij + alpha*sqr(rij)
    //   then s = phi*sqrt(s2)
    double dsds2 = sqr(phi*sigmap)/(2.0*sigma);
    double dsdsp = sigma/sigmap;
    double dsdr = dsds2*(2.0*phase.delta2/cube(dist)
                    + phase.delta1/sqr(dist) + 2.0*sqr(qbroad)*dist);
#else
    //   define s2 = sp**2 - delta2/sqr(rij) - delta1/rij + alpha*sqr(rij)
    //   then s = phi*sqrt(s2)
    double dsds2 = sqr(phi)/(2.0*sigma);
    if (sigma==1e-5) dsds2 = 0;

    double dsdsp = 2.0*sigmap*dsds2;
    double dsdr = dsds2*(2.0*phase.delta2/cube(dist)
                    + phase.delta1/sqr(dist) + 2.0*sqr(qbroad)*dist);
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
        }
        //*pout << phase.a0[i] << " " << d[i] << " " << ig << " " << i << " " << fit_a_i[ipar] << endl;

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
    //------ Derivatives per phase : lat,delta2,delta1
    //------ ----------------------------------------------------------------

    // ----- d/d(lattice parameter a)

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

    // ----- d/d(delta2)


    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = -dTds*dsds2/sqr(dist);
        fit_a_i[ipar] += dg;
    }

    // ----- d/d(delta1)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = -dTds*dsds2/dist;
        fit_a_i[ipar] += dg;
    }

    // ----- d/d(pscale)

    ioffset++;

    // ----- d/d(spdiameter)

    ioffset++;

    // ----- d/d(sratio)

    if ( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = dTds*dsdphi;
        fit_a_i[ipar] += dg;
    }


    //------ ----------------------------------------------------------------
    //------ Derivatives per data set : dscale, qdamp, qbroad
    //------ ----------------------------------------------------------------

    ioffset = this->offset;

    // ----- d/d(dscale[is])

    ioffset++;

    // ----- d/d(qdamp[is])

    ioffset++;

    // ----- d/d(qbroad[ip])

    if( (ipar=fit.refvar[ioffset++]) != -1)
    {
        dg = dTds*dsds2*(2.0*qbroad*sqr(dist));
        fit_a_i[ipar] += dg;
    }

}


/****************************************************************
* Calculated average atom mass in this phase
*****************************************************************/
double Phase::averageAtomicMass()
{
    double mavg = 0.0;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
        mavg += ai->occ * ai->atom_type->M;
    }
    mavg /= np;
    return mavg;
}

/****************************************************************
* Calculated average scattering factor for given radiation type
*****************************************************************/
double Phase::averageScatteringFactor(char tp)
{
    double bavg = 0.0;
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
        bavg += ai->occ * ai->atom_type->sf(tp);
    }
    bavg /= np;
    return bavg;
}

/****************************************************************
* Update atom weights for given scattering type
*****************************************************************/
void Phase::setup_weights(char tp)
{
    // calculate average scattering factor
    double bavg = averageScatteringFactor(tp);
    // get normalized weight of each atom
    for (VAIT ai = atom.begin(); ai != atom.end(); ++ai)
    {
        ai->weight = ai->atom_type->sf(tp) / bavg;
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
    // FIXME - it should be possible to have rcmin smaller than rmin,
    // but there is too much spaghetti that depends on it.
    rcmin = max(rmin, rfmin - rext);
    rcmax = rfmax + rext;
    ncmin = nint((rcmin - rmin)/deltar);
    ncmax = nint((rcmax - rmin)/deltar);
    if (lout)
    {
	*pout << " Extending PDF search distance to " <<
	    rcmin << " -> " << rcmax << " A ...\n";
    }
}

/************************************************************************
 * Diameter of sphere that can enclose primitive cell.
 * This is equal to the longest unit cell diagonal.
 ************************************************************************/
double Phase::circum_diameter()
{
    if (atom.empty())	return 0.0;
    // array of all 4 diagonals
    const size_t numdiags = 4;
    static double ucdiagonals[numdiags][3] = {
        {+1.0, +1.0, +1.0},
        {-1.0, +1.0, +1.0},
        {+1.0, -1.0, +1.0},
        {+1.0, +1.0, -1.0}
    };
    double maxnorm = -1;
    for (size_t idx = 0; idx != numdiags; ++idx)
    {
        const double* ucd = ucdiagonals[idx];
        double normucd = sqrt(skalpro(ucd, ucd));
        if (normucd > maxnorm)
        {
            maxnorm = normucd;
        }
    }
    // adjust to round-off errors
    const double epsilond = sqrt(numeric_limits<double>().epsilon());
    maxnorm = maxnorm*(1.0+epsilond) + epsilond;
    return maxnorm;
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
    }
    ofstream fout;
    stringstream outfilestream;
    string outfilestring = "";

    outfilestream << " " << string(78,'=') << endl
    << " PDF REFINEMENT\n"
    << "   Using PDFFIT version : " << PdfFit::version() << endl
    << " " << string(78,'=') << endl;

    for(int ip=0; ip<nphase; ip++)
        (*phase[ip]).output(outfilestream);

    for (int is=0; is<nset; is++)   datasets[is]->output(outfilestream);

    fit.output(outfilestream);

    outfilestream << " " << string(78,'=') << endl;

    if( fname != "" )
    {
        fout.open(fname.c_str());

        if (!fout)
	{
            throw IOError("Cannot create output file");
        }

        *pout << " Saving fit results to file : " << fname << endl;

        fout << outfilestream.str();
        //fout.setf(ios::showpoint);

        fout.close();
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
        throw unassignedError("data set does not exist");
    }
    else if( fname != "" )
    {
        ofstream fout;
        fout.open(fname.c_str());

        if (!fout)
	{
            throw IOError("cannot create output file");
        }

        *pout << " Saving PDF data set " << iset << " to file : " << fname << endl;

        outfilestring = datasets[iset-1]->build_pdf_file();

        fout << outfilestring;
        fout.close();
    }
    else
    {
        outfilestring = datasets[iset-1]->build_pdf_file();
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
        throw unassignedError("Data set does not exist");
    }
    else if (fname != "" )
    {
        fout.open(fname.c_str());

        if (!fout)
	{
            throw IOError("Cannot create output file");
        }

        *pout << " Saving difference data set " << iset << " to file : " << fname << endl;

        outfilestring = datasets[iset-1]->build_dif_file();

        fout << outfilestring;
        fout.close();
    }
    else
    {
        outfilestring = datasets[iset-1]->build_dif_file();
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
        throw unassignedError("No data set selected");
    }
    assert(nphase == (int)phase.size());
    if (ip == ALL)
    {
        curset->psel = phase;
    }
    else
    {
        // check if one-based index ip is out of bounds
        if (ip < 1 || ip > nphase)
	{
            stringstream eout;
            eout << "Phase " << ip << " undefined";
            throw unassignedError(eout.str());
        }
        curset->selphase(ip - 1, phase[ip - 1]);
    }
}

void DataSet::selphase(int ip, Phase *phase)
{
    if (int(psel.size()) <= ip)
    {
        psel.resize(ip+1);
    }
    psel[ip] = phase;
}


vector<double> DataSet::getcrw() const
{
    assert(mowner);
    double wsqobs = mowner->totalWeighedSquareObs();
    // Get reciprocal value of wsqobs.
    // Do not normalize when wsqobs is zero.
    double recwsqobs = (wsqobs > 0) ? (1.0 / wsqobs) : 1.0;
    vector<double> rv = this->cumchisq;
    vector<double>::iterator xi;
    for (xi = rv.begin(); xi != rv.end(); ++xi)
    {
        *xi = sqrt(*xi * recwsqobs);
    }
    return rv;
}


double DataSet::weighedSquareObs() const
{
    double rv = 0;
    for (int i = nfmin; i <= nfmax; i++)
    {
        rv += wic[i] * obs[i] * obs[i];
    }
    return rv;
}


double DataSet::getdsrw() const
{
    vector<double> crw = this->getcrw();
    double rv = crw.empty() ? 0.0 : crw.back();
    return rv;
}


double DataSet::getdsredchisq() const
{
    assert(mowner);
    int nredobs = mowner->totalReducedObservations();
    double c2 = this->cumchisq.empty() ? 0.0 : this->cumchisq.back();
    double rv = (nredobs > 0) ? (c2 / nredobs) : 0.0;
    return rv;
}


void PdfFit::pdesel(int ip)
{
    if (!curset)
    {
        throw unassignedError("No data set selected");
    }
    assert(nphase == (int)curset->psel.size());
    if (ip == ALL)
    {
        fill(curset->psel.begin(), curset->psel.end(),
                static_cast<Phase*>(NULL));
    }
    else
    {
        // check if one-based index ip is out of bounds
        if (ip < 1 || ip > nphase)
	{
            stringstream eout;
            eout << "phase " << ip << " undefined";
            throw unassignedError(eout.str());
        }
        curset->psel[ip - 1] = NULL;
    }
}

// phase[ip-1] or curphase for ip == 0
Phase* PdfFit::getphase(int ip)
{
    Phase* ph = (0 < ip && ip <= nphase) ? phase[ip-1] : curphase;
    if (!ph || ip < 0 || ip > nphase)
    {
        throw unassignedError("Phase does not exist.");
    }
    return ph;
}

void PdfFit::check_sel_args(int ip, char ijchar, int aidx1)
{
    ostringstream emsg;
    if (!curset)
    {
        throw unassignedError("No data set selected");
    }
    if (ip < 1 || ip > int(curset->psel.size()))
    {
        emsg << "phase " << ip << " undefined or not selected\n";
        throw unassignedError(emsg.str());
    }
    if (ijchar != 'i' && ijchar != 'j')
    {
	ostringstream emsg;
	emsg << "Invalid value of ijchar '" << ijchar << "'";
	throw ValueError(emsg.str());
    }
    if (aidx1 < 1 || aidx1 > (curset->psel[ip - 1]->natoms))
    {
	emsg << "invalid atom index " << aidx1 << ".\n";
	throw ValueError(emsg.str());
    }
}

void PdfFit::selphaseForEachDataSet(Phase* ph)
{
    // find 0-based index of ph in PdfFit::phase vector
    assert(count(this->phase.begin(), this->phase.end(), ph) > 0);
    int phidx0 = find(this->phase.begin(), this->phase.end(), ph) -
        this->phase.begin();
    vector<DataSet*>::iterator dsi;
    for (dsi = this->datasets.begin(); dsi != this->datasets.end(); ++dsi)
    {
        DataSet* pds = *dsi;
        pds->selphase(phidx0, ph);
    }
}

void PdfFit::selectAtomType(int ip, char ijchar, char* symbol, bool select)
{
    check_sel_args(ip, ijchar);
    Phase* ph = curset->psel[ip - 1];
    const LocalPeriodicTable* lpt = ph->getPeriodicTable();
    const AtomType* atp = lpt->lookup(symbol);
    set<int>& ignored = ijchar == 'i' ?
	curset->phase_ignore_i[ph] : curset->phase_ignore_j[ph];
    for (int aidx = 0; aidx < ph->natoms; ++aidx)
    {
	if (atp != ph->atom[aidx].atom_type)	continue;
	if (select)	ignored.erase(aidx);
	else		ignored.insert(aidx);
    }
}

void PdfFit::selectAtomIndex(int ip, char ijchar, int aidx1, bool select)
{
    check_sel_args(ip, ijchar, aidx1);
    Phase* ph = curset->psel[ip - 1];
    set<int>& ignored = ijchar == 'i' ?
	curset->phase_ignore_i[ph] : curset->phase_ignore_j[ph];
    int aidx = aidx1 - 1;
    if (select)	    ignored.erase(aidx);
    else	    ignored.insert(aidx);
}

void PdfFit::selectAll(int ip, char ijchar)
{
    check_sel_args(ip, ijchar);
    Phase* ph = curset->psel[ip - 1];
    set<int>& ignored = ijchar == 'i' ?
	curset->phase_ignore_i[ph] : curset->phase_ignore_j[ph];
    ignored.clear();
}

void PdfFit::selectNone(int ip, char ijchar)
{
    check_sel_args(ip, ijchar);
    Phase* ph = curset->psel[ip - 1];
    set<int>& ignored = ijchar == 'i' ?
	curset->phase_ignore_i[ph] : curset->phase_ignore_j[ph];
    for (int aidx = 0; aidx < ph->natoms; ++aidx)   ignored.insert(aidx);
}

/*****************************************
    Wed Oct 12 2005 - CLF
    Reads observed PDF from arrays.
******************************************/
int PdfFit::read_data_arrays(char tp, double qmax, double qdamp,
        int length, double * r_data, double * Gr_data, double * dGr_data, string _name)
{
    DataSet* pds = new DataSet(this);

    try {
	pds->read_data_arrays(nset+1, tp, qmax, qdamp, length,
		r_data, Gr_data, dGr_data, _name);
    }
    catch(Exception e) {
        delete pds;
        throw;
    }

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        pds->selphase(ip, this->phase[ip]);

    this->datasets.push_back(pds);
    nset++;
    setdata(nset);

    return 1;
}

/*****************************************
    Wed Oct 12 2005 - CLF
    Reads observed PDF from a c-style string.
******************************************/

int PdfFit::read_data_string(string& buffer, char tp, double qmax, double qdamp, string _name)
{
    DataSet* pds = new DataSet(this);
    try {
	pds->read_data_string(nset+1, buffer, tp, qmax, qdamp);
    }
    catch(Exception e) {
	delete pds;
	throw;
    }
    // automatically select existing phases and its atoms for the new dataset
    for (int ip = 0; ip < nphase; ip++)	pds->selphase(ip, this->phase[ip]);

    this->datasets.push_back(pds);
    nset++;
    setdata(nset);

    return 1;
}

/*****************************************
    Reads observed PDF as xy ASCII file.
******************************************/

int PdfFit::read_data(string datafile, char tp, double qmax, double qdamp)
{
    DataSet* pds = new DataSet(this);
    try {
        pds->read_data(nset+1, datafile, tp, qmax, qdamp);
    }
    catch(Exception e) {
        delete pds;
        throw;
    }

    // automatically select existing phases and its atoms for the new dataset
    for (int ip=0; ip<nphase; ip++)
        pds->selphase(ip, this->phase[ip]);

    this->datasets.push_back(pds);
    nset++;
    setdata(nset);

    return 1;
}

// local helper to check for regular spacing in sequence
namespace {

template <class Iterator>
bool isRegular(Iterator first, Iterator last)
{
    if (last - first < 2)    return true;
    double dx = double( *(last-1) - *first ) / double(last - first - 1);
    for (Iterator p0 = first, p1 = first + 1; p1 != last; ++p0, ++p1)
    {
	if (fabs(*p1 - *p0 - dx) > deltar_tol)	    return false;
    }
    return true;
}

}   // local namespace


/* Wed Oct 12 2005 - CLF
 * Using read_data_arrays adds functionality
 * to pdffit2, allowing one to read data that is alread stored as arrays.
 */
void DataSet::read_data_arrays(int _iset, char tp, double _qmax, double _qdamp,
        int length, double * r_data, double * Gr_data,
        double * dGr_data, string _name )
{
    iset = _iset;

    //------ Now analyse given parameters

    //------ - get exp. method (neutron/x-ray)

    scattering_type = tp;

    //------ - get QMAX and Qdamp

    qmax = _qmax;
    qdamp = _qdamp;

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

    *pout << " Reading data from arrays...\n";

    rmin = rfmin = r_data[0];
    rmax = rfmax = r_data[length - 1];
    bin = length;
    deltar = (rmax - rmin)/double(bin-1);
    // check if r has equidistant spacing
    if (!isRegular(r_data, r_data + length))
    {
	throw dataError("Irregular spacing of r values.");
    }
    name = _name;

    *pout << " Read PDF data set " << iset <<
        "  (r = " << rmin << " to " << rmax << " A, " <<
        bin << " points) ...\n";
    if (!lwei)  this->warningOnMissingWeights();
    *pout << endl;

    return;
}

void DataSet::read_data_stream(int _iset, istream& fdata,
	char tp, double _qmax, double _qdamp, string _name)
{
    string line;

    iset = _iset;

    //------ - get exp. method (neutron/x-ray)

    scattering_type = tp;

    //------ - get QMAX and Qdamp

    qmax = _qmax;
    qdamp = _qdamp;

    //------ Ignore header lines
    getline(fdata, line);
    if (line.compare(0, 7, "History") == 0)
    {
	for ( ; !fdata.eof(); getline(fdata,line))
	{
	    if (line[0] != '#')	    continue;
	    // get 3 words from line
	    string w0, w1, w2;
	    istringstream fline(line);
	    fline >> w0 >> w1 >> w2;
	    if (    w0.find_first_not_of('#') == string::npos &&
		    w1 == "start" && w2 == "data" )	    break;
	}
    }
    //------ Any other header lines starting with # ?
    while (line[0] == '#')	getline(fdata,line);

    //------ Finally we actually read the data

    // find number of columns
    int ncol = 0;
    if (!line.empty())
    {
        double x;
        istringstream sline(line);
        while (sline >> x)  ++ncol;
    }

    vector<double> r_data;
    bool lwei = (ncol > 2);     // flag for weights defined by dGr

    while (true)
    {
        double ri, obs;
        double val, wic;
        istringstream sline(line);

        sline >> ri >> obs;
	if (!sline)	break;

        // Obtain weights from dGr.  Use dGr values only when they are all
        // positive, otherwise set all weights to 1.
        wic = 1.0;
        switch (ncol)
        {
            case 3:
                if (sline >> val && val > 0.0)
                {
                    wic = 1.0/sqr(val);
                }
                else
                {
                    lwei = false;
                }
                break;
            case 4:
                // skip one value
                if (sline >> val >> val && val > 0.0)
                {
                    wic = 1.0/sqr(val);
                }
                else
                {
                    lwei = false;
                }
                break;
        }

        // copy values to data arrays
	r_data.push_back(ri);
        this->obs.push_back(obs);
        this->wic.push_back(wic);

        if (!getline(fdata, line))	break;
    }

    // make sure all wic values are one when lwei is false, because
    // lwei could be reset due to zero dGr value
    if (!lwei)
    {
        fill(this->wic.begin(), this->wic.end(), 1.0);
    }

    *pout << " Reading " << ncol << " columns ...\n";

    if (!isRegular(r_data.begin(), r_data.end()))
    {
	throw dataError("Irregular spacing of r values.");
    }
    if (this->obs.size() < 2)
    {
	throw dataError("Incredibly short data set.");
    }
    bin = this->obs.size();
    this->rmin = this->rfmin = r_data.front();
    this->rmax = this->rfmax = r_data.back();
    this->nfmin = 0;
    this->nfmax = bin - 1;
    this->deltar = (rmax - rmin)/double(bin-1);
    this->name = _name;

      *pout << " Read PDF data set " << iset << "  (r = " << rmin
            << " to " << rmax << " A, " << bin << " points) ...\n";

    if (!lwei)  this->warningOnMissingWeights();

    *pout << endl;

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
void DataSet::read_data_string(int _iset, string& buffer, char tp, double _qmax,
        double _qdamp, string _name)
{
    istringstream fdata(buffer);
    read_data_stream(_iset, fdata, tp, _qmax, _qdamp, _name);
    return;
}

/* Wed Oct 12 2005 - CLF
 * Using read_data and the above read_data_string adds functionality
 * to pdffit2, allowing one to read data that has already been loaded.
 */
void DataSet::read_data(int _iset, string pfile, char tp, double _qmax,
        double _qdamp)
{
    // open and check pfile
    ifstream fdata(pfile.c_str());
    if (!fdata)	    throw IOError("File does not exist");
    // read the data
    read_data_stream(_iset, fdata, tp, _qmax, _qdamp, pfile);
    return;
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
	for(is = 0; is < nset; is++)	datasets[is]->range(rmin,rmax);
    }
    else
    {
	if ( (is >= 1) && (is <= nset) )
	{
	    datasets[is-1]->range(rmin,rmax);
	}
	else
	{
	    throw ValueError("Invalid data set number");
	}
    }
}

void DataSet::range(double rmin, double rmax)
{
    if ( (rmin >= this->rmin) && (rmin <= this->rmax)
        && (rmax <= this->rmax) && (rmin < rmax) )
    {
        this->rfmin = rmin;
        this->rfmax = rmax;
    }
    else
    {
        throw ValueError("Range outside data set limits");
    }
}

vector< pair<double,double> >  DataSet::getAtomPhaseFractions()
{
    size_t nphase = psel.size();
    valarray<double> xi(nphase);
    valarray<double> dxi(nphase);
    for (size_t ip = 0; ip < nphase; ip++)
    {
        Phase* ph = psel[ip];
        if (!ph)
        {
            xi[ip] = 0.0;
            dxi[ip] = 0.0;
        }
        else
        {
            double bavg = ph->averageScatteringFactor(scattering_type);
            xi[ip] = ph->pscale / (bavg*bavg);
            dxi[ip] = ph->dpscale / (bavg*bavg);
        }
    }
    double xtot = xi.sum();
    vector< pair<double,double> > rv(nphase, make_pair(0.0, 0.0));
    double dx2tot = (dxi * dxi).sum();
    // get normalized phase fractions fi, do this only when xtot > 0
    for (size_t ip = 0; ip < nphase && 0.0 < xtot; ip++)
    {
        double fi = xi[ip] / xtot;
        double dfi2 = ( dxi[ip]*dxi[ip] * (xtot*xtot - 2*xtot*xi[ip]) +
            xi[ip]*xi[ip]*dx2tot ) / pow(xtot, 4);
        double dfi = sqrt(dfi2);
        rv[ip].first = fi;
        rv[ip].second = dfi;
    }
    return rv;
}

vector< pair<double,double> >  DataSet::getCellPhaseFractions()
{
    size_t nphase = psel.size();
    valarray<double> xi(nphase);
    valarray<double> dxi(nphase);
    for (size_t ip = 0; ip < nphase; ip++)
    {
        Phase* ph = psel[ip];
        if (!ph)
        {
            xi[ip] = 0.0;
            dxi[ip] = 0.0;
        }
        else
        {
            double bavg = ph->averageScatteringFactor(scattering_type);
            xi[ip] = ph->pscale / (bavg*bavg * ph->np);
            dxi[ip] = ph->dpscale / (bavg*bavg * ph->np);
        }
    }
    double xtot = xi.sum();
    vector< pair<double,double> > rv(nphase, make_pair(0.0, 0.0));
    double dx2tot = (dxi * dxi).sum();
    // get normalized phase fractions fi, do this only when xtot > 0
    for (size_t ip = 0; ip < nphase && 0.0 < xtot; ip++)
    {
        double fi = xi[ip] / xtot;
        double dfi2 = ( dxi[ip]*dxi[ip] * (xtot*xtot - 2*xtot*xi[ip]) +
            xi[ip]*xi[ip]*dx2tot ) / pow(xtot, 4);
        double dfi = sqrt(dfi2);
        rv[ip].first = fi;
        rv[ip].second = dfi;
    }
    return rv;
}

vector< pair<double,double> >  DataSet::getMassPhaseFractions()
{
    size_t nphase = psel.size();
    valarray<double> xi(nphase);
    valarray<double> dxi(nphase);
    for (size_t ip = 0; ip < nphase; ip++)
    {
        Phase* ph = psel[ip];
        if (!ph)
        {
            xi[ip] = 0.0;
            dxi[ip] = 0.0;
        }
        else
        {
            double bavg = ph->averageScatteringFactor(scattering_type);
            double mavg = ph->averageAtomicMass();
            xi[ip] = ph->pscale * mavg / (bavg*bavg);
            dxi[ip] = ph->dpscale * mavg / (bavg*bavg);
        }
    }
    double xtot = xi.sum();
    vector< pair<double,double> > rv(nphase, make_pair(0.0, 0.0));
    double dx2tot = (dxi * dxi).sum();
    // get normalized phase fractions fi, do this only when xtot > 0
    for (size_t ip = 0; ip < nphase && 0.0 < xtot; ip++)
    {
        double fi = xi[ip] / xtot;
        double dfi2 = ( dxi[ip]*dxi[ip] * (xtot*xtot - 2*xtot*xi[ip]) +
            xi[ip]*xi[ip]*dx2tot ) / pow(xtot, 4);
        double dfi = sqrt(dfi2);
        rv[ip].first = fi;
        rv[ip].second = dfi;
    }
    return rv;
}


void DataSet::warningOnMissingWeights() const
{
    *pout <<
        " ****WARN****\n" <<
        " Uncertainties on G(r) were absent or unreadable in your input\n" <<
        " data.  The program reset these uncertainties to unity.  This\n" <<
        " does not affect at all the refined parameter values.  However,\n" <<
        " the values of the estimated uncertainties on these refined\n" <<
        " parameter values are not reliable.\n" <<
        " ****WARN****\n";
}

// End of file
