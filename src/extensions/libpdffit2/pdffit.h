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
* classes Phase, DataSet, Fit, Pdf, PdfFit, RefVar, NonRefVar, Builtin
*
* Comments: Main header file included by all others.  Big mess.
*
***********************************************************************/

#ifndef PDFFIT_H_INCLUDED
#define PDFFIT_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "AtomType.h"
#include "Atom.h"
#include "PairDistance.h"
#include "matrix.h"
#include "exceptions.h"
#include "LocalPeriodicTable.h"
#include "OutputStreams.h"

using namespace std;

/***********************************************************************
 *   Here are constants for the parameter coding - DO NOT CHANGE
 *   unless you are rewriting the program !!!!!!!!!!!!!!!!!!!!!!
 ***********************************************************************
 *
 *   n_st    : Number of global structural parameters per phase
 *   n_at    : Number of parameters for each atom
 *   n_ex    : Number of experimental parameters per dataset
 *
 ***********************************************************************/

const int n_st = 11;
const int n_at = 10;
const int n_ex =  3;

const int ALL = -1;
typedef vector<Atom>::iterator VAIT;

enum FCON { USER, IDENT, FCOMP, FSQR };

typedef double (*fbuiltin)(double);

class Fit;
class PdfFit;
class DataSet;
class Phase;
class NonRefVar;
class RefVar;


// non-refinable variables accessible to the users
class NonRefVar
{
    friend class PdfFit;
    double *a;

    public:
    NonRefVar() : a(NULL)
    {
        mtype = "NonRefVar";
    }
    bool isAssigned()  { return (a != NULL); }
    void setptr(double* a)  { this->a = a; }
    void setval(double a)  { *this->a = a; }
    double get() {
        double rv = (a) ? *a : 0.0;
        return rv;
    }
    const string& type() const
    {
        return mtype;
    }

    protected:
    string mtype;
};

// Refinable variables accessible to the users
class RefVar: public NonRefVar
{
    friend class PdfFit;
    public:
    RefVar() : NonRefVar()
    {
        mtype = "RefVar";
    }
};


typedef double (*fcon)(vector<double>&, vector<double>&);
typedef double (*fbuiltin) (double);

// TK commented out 03/31/05
// class Builtin { fbuiltin func, deriv;
//  public:
//  Builtin(fbuiltin func, fbuiltin deriv): func(func), deriv(deriv) {}
// };

// TK 03/31/05 replaced above with
class Builtin {
    public:
    fbuiltin func, deriv;   // made these public.
                            // Really should just be a struct, or provide
                            // const pointers.
    Builtin(fbuiltin func, fbuiltin deriv): func(func), deriv(deriv) {}
};

// contains the address of the variable, and the current value of the parameter
class Fit
{
    friend class PdfFit;

    map<string,Builtin> builtin;

    double parse(string line, vector<double> &dnumdp);
    double compute(string &expr, vector<double> &dnumdp);
    string substitute_pars(string &expr);
    double getnum(istringstream &inexpr, vector<double> &dnumdp);
    void init_builtins();
    void reset();

    public:
    Fit()
    {
	reset();
    }

    // CLF Wed May 25 2005
    // Made these members of fit. Were local variables to
    // refine function.
    double alambda, chisq, ochisq;
    int stagnating;
    // CLF

    int iter;
    double fit_rw;
    double redchisq;  // reduced chi-squared
    double wnorm;   // sum of weighted squared datapoints
    int ntot, ndof;   // total # of points, number of degrees of freedom
    // parameter related variables
    vector<double> p;   // fit-parameters
    vector<double> dp;  // errors on the refined parameters;
    vector<unsigned int> id;     // parameter identification number (arbitrary integer)
    vector<int> ip;     // parameter selection
    matrix<double> covar, alpha;  // covariance matrix and curvature

    // constraint related variables
    vector<double*> var;  // constrained variables
    vector<double> dvar;  // errors on the constrained variables
    vector<bool> vref;    // true if variable contains free parameters
    matrix<double> vcovar;  // covariance matrix on constrained variables
    vector<fcon> fconstraint;  // constraint equations
    vector<string> form;       // constraint formula
    vector<int> idef;  // default parameter if no explicit constraint
    vector<FCON> ctype;   // type of constraint

    vector<int> used;   // vector of used parameter indices (not ids) in the current constraint

    // variables relating constraint-parameter
    matrix<double> dvdp;  // derivative of var wrt p

    // variables relating to refinable variables
    vector<double*> sdptr;  // pointer to standard deviation of refinable variable

    vector<int> refvar;  // integer pointer from refinable variable to actual constraint #

    void fixpar(int pidx);
    void freepar(int pidx);
    void setpar(unsigned int pidx, double val);
    double getpar(unsigned int pidx);
    void constrain(double &a, string form);
    void constrain(double &a, double f(vector<double>&, vector<double>&) );
    void constrain(double &a, int ipar );
    void constrain(double &a, int ipar, FCON type);
    void constrain(double &a, string inpform, fcon f, int idef, FCON type);
    int vfind(double &a);   // look for variable in list of constraints
    int parfind(unsigned int j);
    void fill_variables();
    void fill_errors();
    int varsize() const  { return var.size(); }
    int psize() const  { return p.size(); }
    //Thu Oct 13 2005 - CLF
    void output(ostream &fout);
    void out();

    private:
    vector<int> order_by_id() const;

};

class PdfFit
{
    // class methods

    public:
        static const string& version(const char* ver=NULL);

    private:

	//Struct cr;
	int nphase;
	int total;  // total # of atoms
	vector<Phase*> phase;

	Fit fit;

	// Dataset parameters
	int nset;
	vector<DataSet*> datasets;
	DataSet *curset;

	void init();

    public:

	void fit_setup();
	void fit_errors();
	void fit_theory(bool ldiff, bool lout);
        double totalWeighedSquareObs() const;
        int totalReducedObservations() const;

    private:

	void initarrays();

	void mrqmin(vector<double> &a, vector<int> &ia, matrix<double> &covar,
		matrix<double> &alpha, double &chisq, double &alamda,
		bool deriv);
	void mrqmin(double a[], int ia[], int ma, double **covar,
		double **alpha, double *chisq, double *alamda, bool deriv);
	void mrqcof(double*, int*, int, double**, double*, double*, bool deriv);
	void fill_variables(vector<double> a);
	int getnpar() { return nset*n_ex + nphase*n_st + total*n_at; }

    public:

	Phase* curphase;
	PdfFit();
	~PdfFit();

	void alloc(char tp, double qmax, double qdamp,
		double rmin, double rmax, int bin);
	void calc();
	int read_struct(string fname);  // returns 1:OK, 0:Error
	int read_data(string fname, char tp, double qmax, double qdamp);
	//Wed Oct 12 2005 - CLF
	int read_struct_string(char * buffer);  // returns 1:OK, 0:Error
	int read_data_string(string& buffer, char tp, double qmax, double
		qdamp, string name = "string");
	int read_data_arrays(char tp, double qmax, double qdamp, int length,
		double * r_data, double * Gr_data, double * dGr_data = NULL,
		string name = "array");
	void reset();
	//Thu Oct 13 2005 - CLF
	string save_pdf(int iset, string fname = "");
	string save_dif(int iset, string fname = "");
	string save_res(string fname = "");
	string save_struct(int ip, string strucfile = "");
	string show_struct(int ip);
	//
	int refine(bool deriv, double toler = 0.00000001);
	int refine_step(bool deriv, double toler = 0.00000001);
	double getrw(void)
	{
	    return fit.fit_rw;
	}
	void setpar(unsigned int pidx, double val)
	{
	    fit.setpar(pidx, val);
	}
	void setpar(unsigned int pidx, RefVar v)
	{
	    fit.setpar(pidx, *v.a);
	}
	double getpar(unsigned int pidx)
	{
	    return fit.getpar(pidx);
	}
	void fixpar(int pidx)
	{
	    fit.fixpar(pidx);
	}
	void freepar(int pidx)
	{
	    fit.freepar(pidx);
	}
	void range(int iset, double rmin, double rmax);

	void constrain(RefVar v, double f(vector<double>&, vector<double>&))
	{
	    fit.constrain(*v.a,f);
	}
	void constrain(RefVar v, string form)
	{
	    fit.constrain(*v.a,form);
	}
	void constrain(RefVar v, int ipar)
	{
	    fit.constrain(*v.a,ipar);
	}
	void constrain(RefVar v, int ipar, FCON type)
	{
	    fit.constrain(*v.a,ipar,type);
	}
	void setphase(int ip);
	void setdata(int is);
	void setvar(RefVar v, double a) { v.setval(a); }
	double getvar(RefVar v) { return v.get(); }
	void setvar(NonRefVar v, double a) { v.setval(a); }
	double getvar(NonRefVar v) { return v.get(); }

	void selphase(int ip);
	void pdesel(int ip);
	Phase* getphase(int ip);
        int num_phases()    { return this->nphase; }
        int num_datasets()  { return this->nset; }

    private:

	void check_sel_args(int ip, char ijchar, int aidx1=1);
        void selphaseForEachDataSet(Phase*);

    public:

	void selectAtomType(int ip, char ijchar, char* symbol, bool select);
	void selectAtomIndex(int ip, char ijchar, int aidx1, bool select);
	void selectAll(int ip, char ijchar);
	void selectNone(int ip, char ijchar);

	pair<double,double> bond_angle(int ia, int ja, int ka);
	PairDistance bond_length_atoms(int ia, int ja);
	vector<PairDistance> bond_length_types(string symi, string symj,
		double bmin, double bmax);

	vector<double> getpdf_obs();
	vector<double> getpdf_fit();

	// current phase and set refinable variable pointers
	vector<RefVar> lat, x, y, z,  u11, u22, u33, u12, u13, u23, occ;
	RefVar pscale;
	RefVar spdiameter;
	NonRefVar stepcut;
        RefVar sratio;
        RefVar delta2, delta1;
	RefVar dscale, qdamp, qbroad;
	NonRefVar rcut;
	int getnfmin();
	int getnfmax();
	double getdeltar();
	double getrmin();
	double getrmax();
	vector<double> getcrw() const;
        map<string, vector<double> > getPhaseFractions();
        double get_scat(char tp, string smbpat);
};

class Pdf
{
    public:

	int nfmin, nfmax, ncmin, ncmax;
	double qmax, qdamp, rmin, rmax, deltar;
	double rfmin, rfmax;    // fit range
	double rcmin, rcmax;    // extended calculation range
	double dscale, ddscale, qbroad, dqbroad, dqdamp;

	Pdf()
	{
	    nfmin = nfmax = ncmin = ncmax = 0;
	    qmax = qdamp = rmin = rmax = deltar = 0.0;
	    rfmin = rfmax = rcmin = rcmax = 0.0;
            dscale = 0.0; ddscale = 0.0;
	    qbroad = dqbroad = dqdamp = 0.0;
	}

	vector<double> pdftot;  // total pdf
	matrix<double> calc;  // ?? pdf for each phase and each point in the dataset
};

class DataSet: public Pdf
{

    private:
	int offset;
        const PdfFit* mowner;
	void applyQmaxCutoff(double* y, size_t len);
	void extendCalculationRange(bool lout);
	string selectedAtomsString(int ip, char ijchar);
	void read_data_stream(int iset, istream& fdata,
		char tp, double qmax, double qdamp, string name);
        void warningOnMissingWeights() const;

    public:

	int iset;   // Dataset index
	char scattering_type;
	string name;

	DataSet(const PdfFit* owner) : Pdf(), mowner(owner)
	{
	    dscale = 1.0; ddscale = 0;
            qbroad = dqbroad = 0.0;
	};
	// pdf-related
	void determine(bool ldiff, bool lout, Fit &par);
	void pdf_derivative (Phase& phase,
		const Atom& atomi, const Atom& atomj, double rk, double sigma,
		double sigmap, double dist, double d[3], double ampl,
		double gaus, Fit &fit, double* fit_a_i);

	const vector<double>& getpdf_fit();
	const vector<double>& getpdf_obs();
	//Thu Oct 13 2005 - CLF
	string build_pdf_file();
	string build_dif_file();
	//
	void read_data(int iset, string fname, char tp, double qmax, double qdamp);
	//Wed Oct 12 2005 - CLF
	void read_data_string(int iset, string& buffer, char tp, double qmax, double qdamp,
		string name = "string");
	void read_data_arrays(int iset, char tp, double qmax, double qdamp,
		int length, double * r_data, double * Gr_data, double * dGr_data = NULL, string name = "array");
	//
	//Thu Oct 13 2005 - CLF
	void output(ostream &fout);
	void range(double rmin, double rmax);

	void fit_setup_derivatives(Fit &par);
	void selphase(int ip, Phase *phase);

        vector<double> getcrw() const;
        double weighedSquareObs() const;
        double getdsrw() const;
        double getdsredchisq() const;

        // phase fraction calculations
        vector< pair<double,double> >  getAtomPhaseFractions();
        vector< pair<double,double> >  getCellPhaseFractions();
        vector< pair<double,double> >  getMassPhaseFractions();

	// fit related
	matrix<double> fit_a, fit_b;  // nbin*npar

	int bin;
	vector<double> obs, wic;

	// phase specific information this dataset: selected, allowed atoms
	vector<Phase*> psel;  // phase selection
	// i and j indices to be ignored when calculating PDF
	map<Phase*, set<int> >  phase_ignore_i;
	map<Phase*, set<int> >  phase_ignore_j;
        vector<double> cumchisq;
	friend void PdfFit::fit_setup();

};

class Phase {

    private:

        // data
	string spcgr, name;
	int offset;
        LocalPeriodicTable _local_periodic_table;

	double ar[3], wrez[3], dar[3], dwrez[3];
	double gten[3][3], dgten[3][3];   // tensor and sd
	double rten[3][3], drten[3][3];   // tensor and sd
	double _eps[3][3][3], _reps[3][3][3], _deps[3][3][3], _dreps[3][3][3];
	double &eps(int i, int j, int k) { return _eps[i][j][k]; }
	double &reps(int i, int j, int k) { return _reps[i][j][k]; }
	double &deps(int i, int j, int k) { return _deps[i][j][k]; }
	double &dreps(int i, int j, int k) { return _dreps[i][j][k]; }

        // methods

	set<size_t> selectAtomsOf(string symbol);

	// Added a return value
	string get_scat_string(char tp, const AtomType* atp);

	// shift to equivalent lattice position nearest to the origin
	void make_nearest(double xyz[3]);

    public:

	vector<const AtomType*> atom_types;
	int iphase;
	double cosa, cosb, cosg, sina, sinb, sing;
	double v, dv, vr, dvr;
	int icc[3];
	// Phase has a number of public elements as it is often cross-referenced

	int natoms;   // total number of atoms in multicell icc0 x icc1 x icc2
	int ncatoms;  // number of atoms in one cell

        // IMPORTANT: atom[i].atom_type must be a pointer to AtomType
        // from _local_periodic_table.  Insertion of a new entry to
        // vector atom must be followed by reassign_atom_type(entry);
	vector<Atom> atom;

	double pscale, dpscale;
	double a0[3], win[3], da0[3], dwin[3];
	double np, dnp, rho0, drho0;  // np: total occupance, rho0: number density

	// pdf-related
	double delta2, sratio, rcut;
	double ddelta2, dsratio, delta1, ddelta1;
	double spdiameter, dspdiameter;	// spherical particle diameter
        double stepcut;
	double dnorm, corr_max;


	Phase()
	{
	    pscale = 1.0; dpscale = 0.0;
            sratio=1.0;
	    a0[0] = a0[1] = a0[2] = da0[0] = da0[1] = da0[2] =
		win[0] = win[1] = win[2] = dwin[0] = dwin[1] = dwin[2] =
		delta2 = ddelta2 = dsratio = rcut = 0.0;
	    delta1 = ddelta1 = corr_max = 0.0;
            spdiameter = dspdiameter = 0.0;
            stepcut = 0.0;
	    icc[0] = icc[1] = icc[2] = ncatoms = natoms = 0;
	    spcgr = "P1";
	    name = "UNNAMED";
	}
	inline size_t nscat()
	{
	    return atom_types.size();
	}
	void read_struct(int iphase, string fname);
	void read_struct_string(int iphase, char * buffer);
    private:
	void read_struct_stream(int _iphase, istream& fstruct);
	void read_header(istream &fstruct, bool &ldiscus);
	void read_atoms(istream &fstruct);
        // Reset atom_type to point to an entry from _local_periodic_table
        void reassign_atom_type(Atom* ap);

    public:
	//Thu Oct 13 2005 - CLF
	void output(ostream &fout);
	template <class Stream> void save_struct(Stream &fout);

	void lattice();
	void show_lattice();
	void tensor(double ten[3][3], double vec[3], double win[3]);
	void dtensor(double vec[3], double win[3], double dten[3][3], double dvec[3], double dwin[3]);

	double skalpro(const double h[3], const double k[3]);
	double dskalpro(double h[3], double k[3], double dh[3], double dk[3]);

	double circum_diameter();	// diameter of a sphere enclosing unit cell
	// mean square displacement of 2 atoms
	double msdAtoms(const Atom& ai, const Atom& aj, double* vl);
        // reference to the local periodic table
        LocalPeriodicTable* getPeriodicTable();

	// pdf-related

        double averageAtomicMass();
        double averageScatteringFactor(char tp);
	void setup_weights(char tp);

	pair<double,double> bond_angle(int ia, int ja, int ka);
	PairDistance bond_length_atoms(int ia, int ja);
	vector<PairDistance> bond_length_types(string symi, string symj,
		double bmin, double bmax);
	// Fri Oct 28 2005 - CLF
	// Added a return value
	void show_scat(char tp);
	string get_scat_string(char tp);
	string get_scat_string(char tp, string smbpat);
	void set_scat(char tp, const string& smbpat, double value);
	void reset_scat(const string& smbpat);

	friend class Atom;
	friend class DataSet;
	friend void PdfFit::fit_setup();
	friend void DataSet::fit_setup_derivatives(Fit &par);
	friend void DataSet::determine(bool ldiff, bool lout, Fit &par);
	friend void DataSet::pdf_derivative (Phase& phase,
		const Atom& atomi, const Atom& atomj, double rk, double sigma,
		double sigmap, double dist, double d[3], double ampl,
		double gaus, Fit &fit, double* fit_a_i);
	friend void PdfFit::fit_theory(bool ldiff, bool lout);
};


#endif	// PDFFIT_H_INCLUDED
