// TK: changed "module pdffit" to "module pdffit2"
%module pdffit2
%{
#include "libpdffit2/pdffit.h"
%}
//%include "pdffit.h"


// TK: suggest more descriptive names (sourceType?)
enum Sctp { X, N };
enum FCON { USER, IDENT, FCOMP, FSQR};

// TK: BOO! Preprocessor global is cheesey, even though unused.
#define ALL -1

namespace pdffit
//class PdfFit
{
	//PdfFit();
	//void atom_select(int iset, int lselect);
//	public:
	%rename (pdfrange) range;
	%rename (pfrac) pscale;
	%rename (bang) bond_angle;
	%rename (blen) bond_length;

	void read_struct(char *fname);  // returns 1:OK, 0:Error
	void read_data(char* fname, Sctp t, double qmax, double sigmaq);
	void range(int iset, double rmin, double rmax);
	void reset();
	void alloc(Sctp t, double qmax, double sigmaq, double rmin, double rmax, int bin);
	void calc();
	int refine();
	void save_pdf(int iset, char *fname);
	void save_dif(int iset, char *fname);
	void save_res(char* fname);
	void save_struct(int ip, char *strucfile);
	void show_struct(int ip);
	
	void constrain(RefVar *v, char *form);
	void constrain(RefVar *v, int ipar);
	void constrain(RefVar *v, int ipar, FCON type);

	void setpar(unsigned int n, double val);
	void setpar(unsigned int n, RefVar *v);
	double getpar(unsigned int n);

	void setvar(RefVar *v, double a);
	double getvar(RefVar *v);
	void setvar(NonRefVar *v, double a);
	double getvar(NonRefVar *v);
	
    double getrw(void);
	
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
	double bond_length(int ia, int ja, double bmin, double bmax);
	void show_scat(Sctp type);
	void set_scat(Sctp type, int itype, double len);
	void set_scat(Sctp type, int itype, double a1, double b1, double a2, double b2, 
		double a3, double b3, double a4, double b4, double c);
	void reset_scat(Sctp type, int itype);


	// current phase and set refinable variable pointers
	RefVar *lat(int), *x(int), *y(int), *z(int), *u11(int), *u22(int), *u33(int), 
			*u12(int), *u13(int), *u23(int), *occ(int);
	RefVar *pscale, *srat, *delta, *gamma; 
	RefVar *dscale, *qsig, *qalp;
	NonRefVar *rcut;
	
};

%pythoncode %{
	from math import *
	import signal, os
	import sys
	
	def idle(expr):
		return
	
	
	sys.ps1 = "pdffit2> "
	sys.displayhook = idle
	
	signal.signal(signal.SIGINT,signal.SIG_DFL)
	
	pfrac = _pdffit.cvar.pfrac
	srat = _pdffit.cvar.srat
	delta = _pdffit.cvar.delta
	gamma = _pdffit.cvar.gamma 
	dscale = _pdffit.cvar.dscale
	qsig = _pdffit.cvar.qsig
	qalp = _pdffit.cvar.qalp
	rcut = _pdffit.cvar.rcut
	
%}

