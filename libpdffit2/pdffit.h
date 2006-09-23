#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "matrix.h"
using namespace std;

#define VERSION "2.0 c"
#define DATE "Mar 13 2006"

/*#########################################################################
c   Here are constants for the parameter coding - DO NOT CHANGE
c   unless you are rewriting the program !!!!!!!!!!!!!!!!!!!!!!
c#########################################################################
c
c   n_st    : Number of global structural parameters per phase
c   n_at    : Number of parameters for each atom
c   n_ex    : Number of experimental parameters per dataset
c
c########################################################################*/

const int n_st = 10;
const int n_at = 10;
const int n_ex =  3;

// Wed Nov 16 2005 - CLF
// These macros are only used for debugging and are no longer needed
//#define _p(a) { cout << a << endl; }
//#define _pp(a) { cout << #a << "=" << a << endl; }
//#define _parr(a, n) {  { cout << #a << "=("; for (int i=0; i<n; i++) cout << a[i] << " "; cout << ")" << endl; } }
//#define sgn(x) ((x)>0?1.0:-1.0)

const int ALL = -1;

//#define nint(x) (int(round(x)))
inline int nint(const float x) {
    return (int) rint(x);
}




enum Sctp { X, N };
enum ErrLoc { ER_APPL, ER_COMM };
enum FCON { USER, IDENT, FCOMP, FSQR};

void throw_exception(int errnum, ErrLoc errloc);
void throw_exception(int errnum, string msg);
void warning(string msg);

int strcmp(string str1, string str2, int minlen);
double dget(istringstream &fin);
int iget(istringstream &fin);

// TK added 03/23/05
namespace pdffit
{
    const double pi = 3.1415926535897931;
}
// const double rad=M_PI/180.0;
// const double pi=M_PI, zpi=2.0*M_PI, fpi=4.0*M_PI;

// all these names in the global namespace--not good!
const double rad= pdffit::pi/180.0;
const double pi=pdffit::pi, zpi=2.0*pdffit::pi, fpi=4.0*pdffit::pi;

#include <cmath>
// end TK

const double ln10=log(10.0);

typedef double (*fbuiltin)(double);

double sind(double arg),  cosd(double arg), tand(double arg);
double asind(double arg),  acosd(double arg), atand(double arg);
double exp10(double arg), neg(double);
double dsin(double), dcos(double), dsind(double), dcosd(double), dtan(double),
dtand(double), dasin(double), dacos(double), dasind(double), dacosd(double),
    datan(double), datand(double), dexp(double), dexp10(double), dlog(double),
    dlog10(double), dsqr(double), dcube(double), dsqrt(double), dneg(double);
inline int mod(int i, int j) { return i % j; }
inline double sqr(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }

class Fit;
class PdfFit;
class DataSet;
class Phase;
class NonRefVar;
class RefVar;

class Exception {
    public:
    string msg;
    Exception(string _msg)
        : msg(_msg) {}
    void PrintException() { cout << "Error: " << msg << endl; }
    string GetMsg() { return msg; }
};

//specific exceptions - mimic python names
class ValueError : public Exception
{
    public:
    ValueError(string _msg)
        : Exception( _msg ) {}
};

class unassignedError : public Exception
{
    public:
    unassignedError(string _msg)
        : Exception( _msg ) {}
};

class IOError : public Exception
{
    public:
    IOError(string _msg)
        : Exception( _msg ) {}
};

class structureError : public Exception
{
    public:
    structureError(string _msg)
        : Exception( _msg ) {}
};

class constraintError : public Exception
{
    public:
    constraintError(string _msg)
        : Exception( _msg ) {}
};

class calculationError : public Exception
{
    public:
    calculationError(string _msg)
        : Exception( _msg ) {}
};

class parseError : public Exception
{
    public:
    parseError(string _msg)
        : Exception( _msg ) {}
};

//This one is used internally, and should not make it to the python layer.
class vgetException : public Exception
{
    public:
    vgetException(string _msg)
        : Exception( _msg ) {}
};


string putxdx(double x, double dx);
string cc(double x, double dx);

// non-refinable variables accessible to the users
class NonRefVar
{
    friend class PdfFit;
    double *a;

    public:
    NonRefVar() { a = NULL; }
    bool isAssigned() {
        if(a) return true;
        else return false;
    }
    void setptr(double* a) { this->a = a; }
    void set(double a) { *this->a = a; }
    double get() {
        if (a) {
            return *a;
        }
        else {
            return 0;
        }
    }
};

// Refinable variables accessible to the users
class RefVar: public NonRefVar
{
    friend class PdfFit;
    public:
    RefVar() { NonRefVar(); }
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

    void fixpar(int n);
    void freepar(int n);
    void setpar(unsigned int ip, double val);
    double getpar(unsigned int ip);
    void constrain(double &a, string form);
    void constrain(double &a, double f(vector<double>&, vector<double>&) );
    void constrain(double &a, int ipar );
    void constrain(double &a, int ipar, FCON type);
    void constrain(double &a, string inpform, fcon f, int idef, FCON type);
    int vfind(double &a);   // look for variable in list of constraints
    int parfind(unsigned int j);
    void fill_variables();
    void fill_errors();
    int varsize() { return var.size(); }
    int psize() { return p.size(); }
    //Thu Oct 13 2005 - CLF
    void output(ostream &fout);
    void out();
};

class PdfFit
{
    //Struct cr;
    int nphase;
    int total;  // total # of atoms
    vector<Phase*> phase;

    double xq;  // ??

    bool la;  // ? some crystal propery

    Fit fit;

    // Dataset parameters
    int nset;
    vector<DataSet*> set;
    DataSet *curset;

    void init();
    void fit_reset();

// TK 03/22/05 made these public:
public:
    void fit_setup();
    void fit_errors();
    void fit_theory(bool ldiff, bool lout);
// TK 03/22/05 added private:
private:
    void initarrays();

    void mrqmin(vector<double> &a, vector<int> &ia, matrix<double> &covar,
        matrix<double> &alpha, double &chisq, double &alamda, bool deriv);
    void mrqmin(double a[], int ia[], int ma, double **covar, double **alpha, double *chisq, double *alamda, bool deriv);
    void mrqcof(double*, int*, int, double**, double*, double*, bool deriv);

    void fill_variables(vector<double> a);
    int getnpar() { return nset*n_ex + nphase*n_st + total*n_at; }

    const string version;
    const string cdate;

    public:
    Phase *curphase;
    PdfFit():  version(VERSION), cdate(DATE) {
        init();
        nphase = nset = total = 0;
        la = false; fit.iter = 0;
        curset = NULL; curphase = NULL;
    }
    void alloc(Sctp t, double qmax, double sigmaq, double rmin, double rmax,
            int bin);
    void atom_select(int iset, int lselect);
    void calc();
    int read_struct(string fname);  // returns 1:OK, 0:Error
    int read_data(string fname, Sctp t, double qmax, double sigmaq);
    //Wed Oct 12 2005 - CLF
    int read_struct_string(char * buffer);  // returns 1:OK, 0:Error
    int read_data_string(string& buffer, Sctp t, double qmax, double sigmaq, string name = "string");
    int read_data_arrays(Sctp t, double qmax, double sigmaq,
        int length, double * r_data, double * Gr_data, double * dGr_data = NULL, string name = "array");
    //
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
    void setpar(unsigned int n, double val) { fit.setpar(n, val); }
    void setpar(unsigned int n, RefVar v) { fit.setpar(n, *v.a); }
    double getpar(unsigned int ip) { return fit.getpar(ip); }
    double getrw(void) { return fit.fit_rw; }
    void fixpar(int n) { fit.fixpar(n); }
    void freepar(int n) { fit.freepar(n); }
    void range(int iset, double rmin, double rmax);

    void constrain(RefVar v, double f(vector<double>&, vector<double>&)) { fit.constrain(*v.a,f); }
    void constrain(RefVar v, string form) { fit.constrain(*v.a,form); }
    void constrain(RefVar v, int ipar) { fit.constrain(*v.a,ipar); }
    void constrain(RefVar v, int ipar, FCON type) { fit.constrain(*v.a,ipar,type); }

    void setphase(int ip);
    void setdata(int is);
    void setvar(RefVar v, double a) { v.set(a); }
    double getvar(RefVar v) { return v.get(); }
    void setvar(NonRefVar v, double a) { v.set(a); }
    double getvar(NonRefVar v) { return v.get(); }

    void selphase(int ip);
    void pdesel(int ip);
    void isel(int ip, int i);
    void idesel(int ip, int i);
    void jsel(int ip, int i);
    void jdesel(int ip, int i);

    double bond_angle(int ia, int ja, int ka);
    double bond_length(int ia, int ja);
    void bond_length(int ia, int ja, double bmin, double bmax);

    vector<double> getpdf_obs();
    vector<double> getpdf_fit();

    // current phase and set refinable variable pointers
    vector<RefVar> lat, x, y, z,  u11, u22, u33, u12, u13, u23, occ;
    RefVar pscale, srat, delta, gamma;
    RefVar dscale, qsig, qalp;
    NonRefVar rcut;
    int getnfmin();
    int getnfmax();
    double getdeltar();
    double getrmin();
    double getrmax();
};

class Pdf {
    protected:
    bool lref;    // true if subtract reference pdf
    vector<double> pdfref;  // reference pdf

    public:
    int nfmin, nfmax, ncmin, ncmax;
    double qmax, sigmaq, rmin, rmax, deltar;
    double rfmin, rfmax, rcmin, rcmax, skal;
    double qalp, dqalp, dskal, dsigmaq;
    // Fri Oct 14 2005 - CLF
    // Wrote constructor to initialize all data members.
    Pdf()
    {   nfmin = nfmax = ncmin = ncmax = 0;
        qmax = sigmaq = rmin = rmax = deltar = 0.0;
        rfmin = rfmax = rcmin = rcmax = skal = 0.0;
        qalp = dqalp = dskal = dsigmaq = 0.0;
    }

    vector<double> pdftot;  // total pdf
    matrix<double> calc;  // ?? pdf for each phase and each point in the dataset
    vector<double> sinc;   // jcrb: i think it gets recalculated too often
    vector<double> getpdf_fit() { return pdftot;}

};

class DataSet: public Pdf
{
    int offset;

    // tensor<double> scat; //  (9,0:MAXSCAT,MAXPHA);  now locally stored in setup

    public:
    int iset;   // Dataset index
    bool lxray;
    string name;

    DataSet()
    {   Pdf::Pdf();
        skal=1.0; dskal=0; qalp = dqalp = 0.0; };

    // pdf-related
    void determine(bool ldiff, bool lout, Fit &par);
    void setup_sinc(bool lout);
    void pdf_derivative (Phase &phase,int ia, int ja, double rg, double sigma, double sigmap,
        double dist, double d[3], double ampl,double gaus, Fit &par, double* fit_a_i);

    vector<double> getpdf_fit();
    vector<double> getpdf_obs() {return obs; }
    //Thu Oct 13 2005 - CLF
    string build_pdf_file();
    string build_dif_file();
    //
    void read_data(int iset, string fname, Sctp t, double qmax, double sigmaq, bool lref);
    //Wed Oct 12 2005 - CLF
    void read_data_string(int iset, string& buffer, Sctp t, double qmax, double sigmaq, bool lref,
            string name = "string");
    void read_data_arrays(int iset, Sctp t, double qmax, double sigmaq, bool lref,
            int length, double * r_data, double * Gr_data, double * dGr_data = NULL, string name = "array");
    //
    //Thu Oct 13 2005 - CLF
    void output(ostream &fout);
    void range(double rmin, double rmax);

    void fit_setup_derivatives(Fit &par);
    void selphase(int ip, Phase *phase);
    void selatom(int ip, int i, vector<vector<bool> > &allowed, bool choice);

    // fit related
    matrix<double> fit_a, fit_b;  // nbin*npar

    int bin;
    vector<double> obs, wic;

    // phase specific information this dataset: selected, allowed atoms
    vector<Phase*> psel;  // phase selection
    vector<vector<bool> > allowed_i, allowed_j;   // selected atom type per phase
        // not a matrix because it is not rectangular

    friend void PdfFit::fit_setup();

};


class Atom {
    // as atoms data are often accessed by other functions we make its data public
    int offset;

    public:
    double pos[3], dpos[3];
    double u[6], du[6];
    int iscat;  // atomtype (index in atomlist contained in phase)
    double occ, docc;
    //Wed Oct 12 2005 - CLF
    // Changed input from ifstream to istream
    int read_atom(istream &fstruct, bool ldiscus, void *phase);
    friend void PdfFit::fit_setup();
    friend void DataSet::fit_setup_derivatives(Fit &par);
    friend void DataSet::pdf_derivative (Phase &phase,int ia, int ja, double rg, double sigma, double sigmap,
        double dist, double d[3], double ampl,double gaus, Fit &par, double* fit_a_i);
};

class Phase {
    string spcgr, name;
    int offset;
    vector<vector<double> > Nscatlen;
    vector<vector<double> > Xscatfac;

    double ar[3], wrez[3], dar[3], dwrez[3];
    double gten[3][3], dgten[3][3];   // tensor and sd
    double rten[3][3], drten[3][3];   // tensor and sd
    double _eps[3][3][3], _reps[3][3][3], _deps[3][3][3], _dreps[3][3][3];
    double &eps(int i, int j, int k) { return _eps[i][j][k]; }
    double &reps(int i, int j, int k) { return _reps[i][j][k]; }
    double &deps(int i, int j, int k) { return _deps[i][j][k]; }
    double &dreps(int i, int j, int k) { return _dreps[i][j][k]; }
    vector<double> weight;
    double bave;

    void get_atoms(int type, vector<bool> &latom);
    void get_scat(int i, double* scat_i, bool lxray);
    // Fri Oct 28 2005 - CLF
    // Added a return value
    string show_scat(Sctp type, int itype);

    // shift to equivalent lattice position nearest to the origin
    void make_nearest(double xyz[3]);

    public:
    vector<string> at_lis;
    int iphase;
    double cosa, cosb, cosg, sina, sinb, sing;
    double v, dv, vr, dvr;
    int icc[3];
    // Phase has a number of public elements as it is often cross-referenced
    int nscat;

    int natoms;  // # atoms in structure
    int ncatoms;  // atom/unit cell??
    vector<Atom> atom;
    double skal, dskal;
    double a0[3], win[3], da0[3], dwin[3];
    double np, dnp, rho0, drho0;  // np: total occupance, rho0: number density

    // pdf-related
    double delta, srat, rcut;
    double ddelta, dsrat, gamma, dgamma;
    double dnorm, corr_max;


    Phase() { nscat = 0; skal=1.0; srat=1.0;
          dskal = a0[1] = a0[1] = a0[2] = da0[0] = da0[1] = da0[2] =
          win[0] = win[1] = win[2] = dwin[0] = dwin[1] = dwin[2] =
          //ar[0] = ar[1] = ar[2] = dar[0] = dar[1] = dar[2] = wrez[0] = wrez[1] =
          //wrez[2] = dwrez[0] = dwrez[1] = dwrez[2] =
          delta = ddelta = dsrat = rcut = 0.0;  // pdf-initializations
          gamma = dgamma = corr_max = 0.0;
          icc[0] = icc[1] = icc[2] = ncatoms = natoms = 0;
          spcgr = "P1";
          name = "UNNAMED";  }
    int getnscat() { return nscat; }
    void setatom(string type) { at_lis.push_back(type); nscat++; }
    void printatoms() { int i; for (i=0; i<nscat; i++) { cout << i+1 << " " << at_lis[i] << endl; } }
    void add_atom();
    int get_iscat(string type);
    void list_atom_types() { int i; for (i=0; i<nscat; i++) cout << at_lis[i] << endl; };

    void read_struct(int iphase, string fname, int &total);
    //Wed Oct 12 2005 - CLF
    void read_struct_string(int iphase, char * buffer, int &total);
    // Changed input from ifstream to istream
    void read_header(istream &fstruct, bool &ldiscus);
    void read_atoms(istream &fstruct, bool ldiscus);
    //
    //Thu Oct 13 2005 - CLF
    void output(ostream &fout);
    template <class Stream> void save_struct(Stream &fout);

    void lattice(bool lout);
    void tensor(double ten[3][3], double vec[3], double win[3]);
    void dtensor(double vec[3], double win[3], double dten[3][3], double dvec[3], double dwin[3]);

    double skalpro(double h[3], double k[3]);
    double dskalpro(double h[3] ,double k[3], double dh[3], double dk[3]);

    double circum_diameter();	// diameter of a sphere enclosing unit cell
    // mean square displacement of 2 atoms
    double msdAtoms(const Atom& ai, const Atom& aj, double* vl);

    // pdf-related

    void setup_weights(bool lout, bool lxray);
    void dlink(matrix<double> &scat, bool lxray);
    double scatteringFactor(double* scat_i, bool lxray);

    double bond_angle(int ia, int ja, int ka);
    double bond_length(int ia, int ja);
    void bond_length(int ia, int ja, double bmin, double bmax);
    // Fri Oct 28 2005 - CLF
    // Added a return value
    string show_scat(Sctp type);
    void set_scat(Sctp type, int itype, double len);
    void set_scat(Sctp type, int itype, double a1, double b1, double a2, double b2,
    double a3, double b3, double a4, double b4, double c);
    void reset_scat(Sctp type, int itype);

    friend class Atom;
    friend class DataSet;
    friend void PdfFit::fit_setup();
    friend void DataSet::fit_setup_derivatives(Fit &par);
    friend void DataSet::determine(bool ldiff, bool lout, Fit &par);
    friend void DataSet::pdf_derivative (Phase &phase,int ia, int ja, double rg, double sigma, double sigmap,
        double dist, double d[3], double ampl,double gaus, Fit &par, double* fit_a_i);
    friend void PdfFit::fit_theory(bool ldiff, bool lout);
};


namespace pdffit
{
    void read_struct(char *fname);  // returns 1:OK, 0:Error
    void read_data(char* fname, Sctp t, double qmax, double sigmaq);
    void range(int iset, double rmin, double rmax);
    void reset();
    void alloc(Sctp t, double qmax, double sigmaq, double rmin, double rmax, int bin);
    void calc();
    int refine();
    int refine(bool deriv, double toler = 0.00000001);
    //Thu Oct 13 2005 - CLF
    string save_pdf(int iset, string fname);
    string save_dif(int iset, string fname);
    string save_res(string fname);
    string save_struct(int ip, string strucfile);
    string show_struct(int ip);
    //

    void constrain(RefVar *v, char *form);
    void constrain(RefVar *v, int ipar);
    void constrain(RefVar *v, int ipar, FCON type);

    void setpar(unsigned int n, double val);
    void setpar(unsigned int n, RefVar *v);
    double getpar(unsigned int ip);
    double getrw(void);

    void setvar(RefVar *v, double a);
    double getvar(RefVar *v);
    void setvar(NonRefVar *v, double a);
    double getvar(NonRefVar *v);

    void fixpar(int n);
    void freepar(int n);

    void setphase(int ip);
    void setdata(int is);

    void psel(int ip);
    void pdesel(int ip);
    void isel(int ip, int i);
    void idesel(int ip, int i);
    void jsel(int ip, int i);
    void jdesel(int ip, int i);

    double bond_angle(int ia, int ja, int ka);
    double bond_length(int ia, int ja);
    void bond_length(int ia, int ja, double bmin, double bmax);
    // Fri Oct 28 2005 - CLF
    // Added a return value
    string show_scat(Sctp type);
    void set_scat(Sctp type, int itype, double len);
    void set_scat(Sctp type, int itype, double a1, double b1, double a2, double b2,
        double a3, double b3, double a4, double b4, double c);
    void reset_scat(Sctp type, int itype);


    // current phase and set refinable variable pointers
    RefVar *lat(int), *x(int), *y(int), *z(int), *u11(int), *u22(int), *u33(int),
                *u12(int), *u13(int), *u23(int), *occ(int);
    extern RefVar *pscale, *srat, *delta, *gamma;
    extern RefVar *dscale, *qsig, *qalp;
    extern NonRefVar *rcut;
}
