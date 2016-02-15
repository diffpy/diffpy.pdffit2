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
* PdfFit methods for Levenberg-Marquardt least-squares minimization.
*
* Comments:
*
***********************************************************************/

//#define CHECK_DERIVATIVES

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "MathUtils.h"
#include "matrix.h"
#include "nrutil.h"
#include "pdffit.h"

using NS_PDFFIT2::pout;

void PdfFit::mrqmin(vector<double> &a, vector<int> &ia, matrix<double> &covar,
	matrix<double> &alpha, double &chisq, double &alamda, bool deriv)
{
    using std::copy;

    int ma = a.size();
    double** _covar = new double*[ma];
    double** _alpha = new double*[ma];
    double* _a = new double[ma];

    copy(a.begin(), a.end(), _a);

    for (int i = 0; i < ma; i++)
    {
	_covar[i] = &covar[i][0]-1;
	_alpha[i] = &alpha[i][0]-1;
    }

    int* iaptr = ia.size() ? &(ia[0]) - 1 : NULL;
    mrqmin(_a-1, iaptr, ma, _covar-1, _alpha-1, &chisq, &alamda, deriv);

    copy(_a, _a + ma, a.begin());

    delete [] _a;
    delete [] _alpha;
    delete [] _covar;
}


// Levenberg-Marquardt method, attempting to reduce the value Chi2 of a fit
// between a set of data points x[1..ndata], y[1..ndata] with individual
// standard deviations sig[1..ndata], and a nonlinear function dependent on ma
// coeffcients a[1..ma].  The input array ia[1..ma] indicates by nonzero
// entries those components of a that should be fitted for, and by zero entries
// those components that should be held fixed at their input values. The
// program re- turns current best-fit values for the parameters a[1..ma], and
// Chi2=chisq. The arrays covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as
// working space during most iterations. Supply a routine
// funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, and its
// derivatives dyda[1..ma] with respect to the fitting parameters a at x. On
// the first call provide an initial guess for the parameters a, and set
// alamda<0 for initialization (which then sets alamda=.001). If a step
// succeeds chisq becomes smaller and alamda de- creases by a factor of 10. If
// a step fails alamda grows by a factor of 10. You must call this routine
// repeatedly until convergence is achieved. Then, make one final call with
// alamda=0, so that covar[1..ma][1..ma] returns the covariance matrix, and
// alpha the curvature matrix.  (Parameters held fixed will return zero
// covariances.)
void PdfFit::mrqmin(double a[], int ia[], int ma, double **covar, double **alpha, double *chisq, double *alamda, bool deriv)
{
    void covsrt(double **covar, int ma, int ia[], int mfit);
    void gaussj(double **a, int n, double **b, int m);

    int j, k, l;
    static int mfit;
    static double ochisq, *atry, *beta, *da, **oneda;

    if (*alamda < 0.0)
    {
        // Initialization.
	atry = dvector(1, ma);
	beta = dvector(1, ma);
	da = dvector(1, ma);

	for (mfit = 0, j = 1; j <= ma; j++)
	    if (ia[j]) mfit++;

	oneda = dmatrix(1, mfit, 1, 1);
	*alamda = 0.001;
	mrqcof(a, ia, ma, alpha, beta, chisq, deriv);
	ochisq=(*chisq);
	for (j = 1; j <= ma; j++) atry[j]=a[j];

	*pout << "\n******************************** ITER: " << fit.iter << " ********************************\n";
	fit.fit_rw = sqrt(ochisq/fit.wnorm);
	fit.redchisq = ochisq/(fit.ntot-fit.ndof);
	fit.out();
	*pout << " chisq.: " << ochisq << "   red.chisq.: " << fit.redchisq << "   Rw: " << fit.fit_rw << endl;
    }
    // Alter linearized fitting matrix, by augmenting diagonal elements.
    for (j = 1; j <= mfit; j++)
    {
	for (k = 1; k <= mfit; k++) covar[j][k]=alpha[j][k];
	covar[j][j]=alpha[j][j]*(1.0+(*alamda));
	oneda[j][1]=beta[j];
    }

    //===================================================================
#if defined(TEST)
    print(covar, mfit); double **save = dmatrix(1, mfit, 1, mfit);
    for(int i = 1; i <= mfit; i++) for(j = 1; j <= mfit; j++) save[i][j] = covar[i][j];
    for(int i = 1; i <= mfit; i++) printf("%lg ", oneda[i][1]); printf("\n");
#endif
    //===================================================================

    gaussj(covar, mfit, oneda, 1); // Matrix solution.

    //===================================================================
#if defined(TEST)
    double res;
    for(int i = 1; i <= mfit; i++)
    {
	res = 0;
	for(j = 1; j <= mfit; j++)
	    res += save[i][j]*oneda[j][1]; printf("%lg ", res);
    }
    printf("\n");
#endif
    //=================================================================

    for (j = 1; j <= mfit; j++) da[j]=oneda[j][1];
    if (*alamda == 0.0) { // Once converged, evaluate covariance matrix.
	covsrt(covar, ma, ia, mfit);
	covsrt(alpha, ma, ia, mfit); // Spread out alpha to its full size too.
	free_dmatrix(oneda, 1, mfit, 1, 1);
	free_dvector(da, 1, ma);
	free_dvector(beta, 1, ma);
	free_dvector(atry, 1, ma);
	return;
    }
    for (j = 0, l = 1; l <= ma; l++) // Did the trial succeed?
	if (ia[l]) atry[l]=a[l]+da[++j];
    mrqcof(atry, ia, ma, covar, da, chisq, deriv);
    if (*chisq < ochisq) { // Success, accept the new solution.
	*alamda *= 0.1;
	ochisq=(*chisq);
	for (j = 1; j <= mfit; j++) {
	    for (k = 1; k <= mfit; k++) alpha[j][k]=covar[j][k];
	    beta[j]=da[j];
	}
	for (l = 1; l <= ma; l++) a[l]=atry[l];
    } else { // Failure, increase alamda and return.
	*alamda *= 10.0;
	*chisq = ochisq;
    }
}

// Used by mrqmin to evaluate the linearized fitting matrix alpha, and
// vector beta as in (15.5.8), and calculate Chi2.
void PdfFit::mrqcof(double a[], int ia[], int ma, double **alpha, double beta[], double *chisq, bool deriv)
{
    int i, j, k, l, m, mfit = 0;
    double wt, sig2i, dy, *dyda;
    //double ymod;

    dyda = dvector(1, ma);
    for (j = 1; j <= ma; j++)
	if (ia[j]) mfit++;
    for (j = 1; j <= mfit; j++) { // Initialize (symmetric) alpha, beta.
	for (k = 1; k <= j; k++) alpha[j][k]=0.0;
	beta[j]=0.0;
    }

    // careful: a of mrqcof is in fact atry of mrqmin!
    for (j = 1; j <= ma; j++)
	fit.p[j-1] = a[j];

    // printf("a=%.12lg %.12lg %.12lg\n", a[1], a[2], a[3]);

    fit_theory(deriv, false);  // yields pdftot and derivatives wrt parameters

    //=============================================================================
#if !defined(CHECK_DERIVATIVES)
    // deriv: true for analytic derivatives, false for numerical ones
    if (!deriv)
#endif
    {
#if defined(CHECK_DERIVATIVES)
	matrix<double> dersave = datasets[0]->fit_b;
#endif

	// numerical derivative
	vector<vector<double> > pdfsave;
	for (int is = 0; is < nset; is++)
	{
	    pdfsave.push_back(datasets[is]->pdftot);

	    // fit_b used both in numerical and analytical derivative case
	    datasets[is]->fit_b.clear();
	    datasets[is]->fit_b.resize(datasets[is]->ncmax+1, fit.psize());
	}

	double psave;

	if (fit.psize() != ma) {
	    throw constraintError("Inconsistent number of parameters");
	}

	for (int ip = 0; ip < fit.psize(); ip++)
	{
	    if (!fit.ip[ip]) continue;

	    double delta = double_eps*fit.p[ip];
	    if (fabs(delta) < double_eps/100)	delta = double_eps;
	    psave = fit.p[ip];
	    fit.p[ip] += delta;
	    fit_theory(false, false);  // yields pdftot and derivatives wrt parameters
	    fit.p[ip] = psave;

	    for (int is = 0; is < nset; is++)
	    {
		DataSet* pds = this->datasets[is];

		for (i = pds->nfmin; i <= pds->nfmax; i++)
		{
		    pds->fit_b[i][ip] = (pds->pdftot[i]-pdfsave[is][i])/delta;
		}
	    }

#if defined(CHECK_DERIVATIVES)
	    i = 200;
	    {
		*pout << "DERIVATIVES:: ANALYTIC : " << dersave[i][ip] << endl;
		*pout << "              NUMERICAL: " << (datasets[0]->pdftot[i]-pdfsave[0][i])/delta
		    << " (delta[" << fit.id[ip] << "]=" << delta << ")" << endl << endl;
	    }
#endif
	}
	for (int is = 0; is < nset; is++)	datasets[is]->pdftot = pdfsave[is];

#if defined(CHECK_DERIVATIVES)
	exit(0);
#endif
    }
    //=============================================================================

    *chisq = 0.0;

    for (int is = 0; is < nset; is++)
    {
	DataSet* pds = this->datasets[is];
        pds->cumchisq.clear();

	for (i = pds->nfmin; i <= pds->nfmax; i++)
	{ // Summation loop over all data.

	    //(*funcs)(x[i], a, &ymod, dyda, ma);
	    sig2i = pds->wic[i];
	    dy = pds->obs[i] - pds->pdftot[i];

	    for (j = 1; j <= ma; j++)
		dyda[j] = pds->fit_b[i][j-1];   // of course use fit_b, NOT fit_a

	    //if (i==0) printf("dyda=%.12lg %.12lg %.12lg\n", dyda[1], dyda[2], dyda[3]);

	    for (j = 0, l = 1; l <= ma; l++) {
		if (ia[l]) {
		    wt = dyda[l]*sig2i;
		    for (j++, k = 0, m = 1; m <= l; m++)
			if (ia[m]) alpha[j][++k] += wt*dyda[m];
		    beta[j] += dy*wt;
		}
	    }
            double prev = pds->cumchisq.empty() ? 0.0 : pds->cumchisq.back();
            double chisq_contribution = dy * dy * sig2i;
            pds->cumchisq.push_back(prev + chisq_contribution);
	}
        *chisq += pds->cumchisq.empty() ? 0.0 : pds->cumchisq.back();
    }
    for (j = 2; j <= mfit; j++) // Fill in the symmetric side.
	for (k = 1; k < j; k++) alpha[k][j]=alpha[j][k];
    free_dvector(dyda, 1, ma);
    //print(alpha, ma);
}


// Expand in storage the covariance matrix covar, so as to take into account
// parameters that are being held fixed. (For the latter, return zero
// covariances.)
void covsrt(double **covar, int ma, int ia[], int mfit)
{
    int i, j, k;
    for (i = mfit+1; i <= ma; i++)
    {
	for (j = 1; j <= i; j++)
	{
	    covar[i][j]=covar[j][i]=0.0;
	}
    }
    k = mfit;
    for (j = ma; j >= 1; j--) {
	if (ia[j])
	{
	    for (i = 1; i <= ma; i++)	swap(covar[i][k], covar[i][j]);
	    for (i = 1; i <= ma; i++)	swap(covar[k][i], covar[j][i]);
	    k--;
	}
    }
}

// End of file
