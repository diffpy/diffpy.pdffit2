/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Bindings from python to c++ PdfFit class.
*
* Comments:
*
***********************************************************************/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "misc.h"
#include "pyexceptions.h"
#include "PyFileStreambuf.h"
#include "../libpdffit2/StringUtils.h"
#include "../libpdffit2/LocalPeriodicTable.h"
#include "../libpdffit2/pdffit.h"

// ostream buffer used for engine output redirection
PyFileStreambuf* py_stdout_streambuf = NULL;

// copyright
char pypdffit2_copyright__doc__[] = "";
char pypdffit2_copyright__name__[] = "copyright";

static char pypdffit2_copyright_note[] =
    "pdffit2 python module: Copyright (c) 2005-2016 Simon J. L. Billinge et al.";

// constant strings for python capsule names (cn)
static char* cnpfit = "pdffit";
static char* cnvar = "pdfvar";

PyObject * pypdffit2_copyright(PyObject *, PyObject *)
{
    return Py_BuildValue("s", pypdffit2_copyright_note);
}

//helper function to convert a pylist to a double array.
void double_array_from_pylist(PyObject *pylist, double *d_array, int const length)
{
    //length is the size of the d_array and not necessarily equal to the length
    //of the pylist
    PyObject *pyval = 0;

    for(int i = 0; i < length; i++) {
        pyval = PyList_GetItem(pylist, i);
        d_array[i] = PyFloat_AsDouble(pyval);
    }
}

// helper function to delete PdfFit object
static void deletePdfFit(PyObject* ptr)
{
    PdfFit *pdf = (PdfFit *)PyCapsule_GetPointer(ptr, cnpfit);
    delete pdf;
    return;
}

// create a PdfFit instance
char pypdffit2_create__doc__[] = "";
char pypdffit2_create__name__[] = "create";

PyObject * pypdffit2_create(PyObject *, PyObject *args)
{
    PdfFit *ppdf = new PdfFit();
    PyObject *py_ppdf = PyCapsule_New((void *)ppdf, cnpfit, deletePdfFit);
    return py_ppdf;
}


// read_struct
char pypdffit2_read_struct__doc__[] = "Read structure file into memory.";
char pypdffit2_read_struct__name__[] = "read_struct";

PyObject * pypdffit2_read_struct(PyObject *, PyObject *args)
{
    char *fname;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Os", &py_ppdf, &fname);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->read_struct(fname);
    }
    catch(structureError e) {
        PyErr_SetString(pypdffit2_structureError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(pypdffit2_structureError, e.GetMsg().c_str());
        return 0;
    }
    catch(calculationError e) {
        PyErr_SetString(pypdffit2_calculationError, e.GetMsg().c_str());
        return 0;
    }
    catch(IOError e) {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// read_struct_string
char pypdffit2_read_struct_string__doc__[] = "Read structure file into memory from a c-string.";
char pypdffit2_read_struct_string__name__[] = "read_struct_string";

PyObject * pypdffit2_read_struct_string(PyObject *, PyObject *args)
{
    char *buffer;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Os", &py_ppdf, &buffer);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->read_struct_string(buffer);
    }
    catch(structureError e) {
        PyErr_SetString(pypdffit2_structureError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(pypdffit2_structureError, e.GetMsg().c_str());
        return 0;
    }
    catch(calculationError e) {
        PyErr_SetString(pypdffit2_calculationError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}


// read_data
char pypdffit2_read_data__doc__[] = "Read data file into memory.";
char pypdffit2_read_data__name__[] = "read_data";

PyObject * pypdffit2_read_data(PyObject *, PyObject *args)
{
    char *fname;
    char stype;
    double qmax, qdamp;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oscdd", &py_ppdf, &fname, &stype, &qmax, &qdamp);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->read_data(fname, stype, qmax, qdamp);
    }
    catch(IOError e) {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(dataError e) {
        PyErr_SetString(pypdffit2_dataError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// read_data_string
char pypdffit2_read_data_string__doc__[] = "Read data from string into memory.";
char pypdffit2_read_data_string__name__[] = "read_data_string";

PyObject * pypdffit2_read_data_string(PyObject *, PyObject *args)
{
    char *buffer;
    char *c_name = NULL;
    char stype;
    double qmax, qdamp;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oscdd|s", &py_ppdf, &buffer, &stype, &qmax, &qdamp, &c_name);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    string name = c_name ? c_name : "";
    try {
	string sbuffer(buffer);
        ppdf->read_data_string(sbuffer, stype, qmax, qdamp, name);
    }
    catch(IOError e) {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(dataError e) {
        PyErr_SetString(pypdffit2_dataError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// read_data_arrays - read_data_lists in PdfFit class
char pypdffit2_read_data_arrays__doc__[] = "Read data from arrays into memory.";
char pypdffit2_read_data_arrays__name__[] = "read_data_arrays";

PyObject * pypdffit2_read_data_arrays(PyObject *, PyObject *args)
{
    char stype;
    double qmax, qdamp;
    int length;
    char * c_name = NULL;
    double *r_data = NULL;
    double *Gr_data = NULL;
    double *dGr_data = NULL;
    PyObject *py_r_data = Py_None;
    PyObject *py_Gr_data = Py_None;
    PyObject *py_dGr_data = Py_None;
    PyObject *py_ppdf = NULL;
    int ok = PyArg_ParseTuple(args, "OcddOO|Os", &py_ppdf, &stype, &qmax, &qdamp,
            &py_r_data, &py_Gr_data, &py_dGr_data, &c_name);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);

    length = PyList_Size(py_Gr_data);
    //quick check that the arrays are all the same length //
    int r_len = PyList_Size(py_r_data);
    int dGr_len = length;
    if(py_dGr_data != Py_None) {
        dGr_len = PyList_Size(py_dGr_data);
    }
    if(r_len != length || dGr_len != length)
    {
	string err_string = "Data arrays have different lengths";
	PyErr_SetString(PyExc_ValueError, err_string.c_str());
	return 0;
    }

    // create data arrays
    r_data = new double [length];
    double_array_from_pylist(py_r_data, r_data, length);
    Gr_data = new double [length];
    double_array_from_pylist(py_Gr_data, Gr_data, length);
    if(py_dGr_data != Py_None) {
        dGr_data = new double [length];
        double_array_from_pylist(py_dGr_data, dGr_data, length);
    }
    string name = c_name;
    try {
	ppdf->read_data_arrays(stype, qmax, qdamp, length,
		r_data, Gr_data, dGr_data, name);
    }
    catch(dataError e) {
        PyErr_SetString(pypdffit2_dataError, e.GetMsg().c_str());
        return 0;
    }

    // read_data_arrays creates its own copy of the data, so we must delete our
    // copies here.
    delete [] r_data;
    delete [] Gr_data;
    if( dGr_data != NULL ) delete [] dGr_data;
    Py_INCREF(Py_None);
    return Py_None;
}

// pdfrange (range in c)
char pypdffit2_pdfrange__doc__[] = "Set r-range of pdf.";
char pypdffit2_pdfrange__name__[] = "pdfrange";

PyObject * pypdffit2_pdfrange(PyObject *, PyObject *args)
{
    int iset;
    double rmin, rmax;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oidd", &py_ppdf, &iset, &rmin, &rmax);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->range(iset, rmin, rmax);
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// reset
char pypdffit2_reset__doc__[] = "reset pdf data";
char pypdffit2_reset__name__[] = "reset";

PyObject * pypdffit2_reset(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    ppdf->reset();
    Py_INCREF(Py_None);
    return Py_None;
}

// alloc
char pypdffit2_alloc__doc__[] = "Allocate space for pdf data";
char pypdffit2_alloc__name__[] = "alloc";

PyObject * pypdffit2_alloc(PyObject *, PyObject *args)
{
    char stype;
    double qmax, qdamp, rmin, rmax;
    int bin;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ocddddi", &py_ppdf, &stype, &qmax, &qdamp, &rmin, &rmax, &bin);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->alloc(stype, qmax, qdamp, rmin, rmax, bin);
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    catch (unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// calc
char pypdffit2_calc__doc__[] = "calculate pdf from data";
char pypdffit2_calc__name__[] = "calc";

PyObject * pypdffit2_calc(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->calc();
    }
    catch (calculationError e) {
        PyErr_SetString(pypdffit2_calculationError, e.GetMsg().c_str());
        return 0;
    }
    catch (unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch (parseError e) {
        PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// refine
char pypdffit2_refine__doc__[] = "refine model to pdf data";
char pypdffit2_refine__name__[] = "refine";

PyObject * pypdffit2_refine(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    double toler;
    int ok = PyArg_ParseTuple(args, "Od", &py_ppdf, &toler);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        ppdf->refine(true, toler);
    }
    catch(constraintError e) {
        PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
        return 0;
    }
    catch(calculationError e) {
        PyErr_SetString(pypdffit2_calculationError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// local helper class for background thread in pypdffit2_refine_step()
namespace {

class RefineStepHelper
{
    private:
	PyThreadState* thread_state;
	ostringstream  msgout;

    public:
	// Constructor saves thread state and arranges for holding
	// engine output when redirected
	RefineStepHelper()
	{
	    if (py_stdout_streambuf)
	    {
		NS_PDFFIT2::pout->rdbuf(msgout.rdbuf());
	    }
	    thread_state = PyEval_SaveThread();
	}

	// method for restoring thread state and writing any outstanding output
	void clean()
	{
	    PyEval_RestoreThread(thread_state);
	    if (py_stdout_streambuf)
	    {
		NS_PDFFIT2::pout->rdbuf(py_stdout_streambuf);
		*NS_PDFFIT2::pout << msgout.str();
	    }
	}
};

}   // local namespace

// refine_step
char pypdffit2_refine_step__doc__[] = "Make one step in the refinement process.";
char pypdffit2_refine_step__name__[] = "refine_step";

PyObject * pypdffit2_refine_step(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    double toler;
    int ok = PyArg_ParseTuple(args, "Od", &py_ppdf, &toler);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    RefineStepHelper janitor;   // takes care of thread an output issues
    int finished = 1;
    try {
	finished = ppdf->refine_step(true, toler);
    }
    catch(parseError e) {
	// parseError is due to invalid constraint
	janitor.clean();
	PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
	return 0;
    }
    catch(constraintError e) {
	janitor.clean();
	PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
	return 0;
    }
    catch(calculationError e) {
	janitor.clean();
	PyErr_SetString(pypdffit2_calculationError, e.GetMsg().c_str());
	return 0;
    }
    catch(unassignedError e) {
	janitor.clean();
	PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
	return 0;
    }
    catch(...) {
	janitor.clean();
	return 0;
    }
    janitor.clean();

    return Py_BuildValue("i", finished);
}

// save_pdf
char pypdffit2_save_pdf__doc__[] = "Save calculated pdf to file";
char pypdffit2_save_pdf__name__[] = "save_pdf";

PyObject * pypdffit2_save_pdf(PyObject *, PyObject *args)
{
    char *fname;
    int iset;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ois", &py_ppdf, &iset, &fname);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        string outfilestring = ppdf->save_pdf(iset, fname);
        return Py_BuildValue("s", outfilestring.c_str());
    }
    catch(IOError e)
    {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// save_dif
char pypdffit2_save_dif__doc__[] = "Save pdf difference to file";
char pypdffit2_save_dif__name__[] = "save_dif";

PyObject * pypdffit2_save_dif(PyObject *, PyObject *args)
{
    char *fname;
    int iset;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ois", &py_ppdf, &iset, &fname);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        string outfilestring = ppdf->save_dif(iset, fname);
        return Py_BuildValue("s", outfilestring.c_str());
    }
    catch(IOError e)
    {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// save_res
char pypdffit2_save_res__doc__[] = "Save residual to file";
char pypdffit2_save_res__name__[] = "save_res";

PyObject * pypdffit2_save_res(PyObject *, PyObject *args)
{
    char *fname;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Os", &py_ppdf, &fname);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        string outfilestring = ppdf->save_res(fname);
        return Py_BuildValue("s", outfilestring.c_str());
    }
    catch(IOError e)
    {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// save_struct
char pypdffit2_save_struct__doc__[] = "Save refined structure to file";
char pypdffit2_save_struct__name__[] = "save_struct";

PyObject * pypdffit2_save_struct(PyObject *, PyObject *args)
{
    char *fname;
    int iset;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ois", &py_ppdf, &iset, &fname);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        string outfilestring = ppdf->save_struct(iset, fname);
        return Py_BuildValue("s", outfilestring.c_str());
    }
    catch(IOError e)
    {
        PyErr_SetString(PyExc_IOError, e.GetMsg().c_str());
        return 0;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// show_struct
char pypdffit2_show_struct__doc__[] = "Show structure.";
char pypdffit2_show_struct__name__[] = "show_struct";

PyObject * pypdffit2_show_struct(PyObject *, PyObject *args)
{
    int ip;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        string outfilestring = ppdf->show_struct(ip);
        return Py_BuildValue("s", outfilestring.c_str());
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// constrain to string
char pypdffit2_constrain_str__doc__[] = "Constrain refined variable to string.";
char pypdffit2_constrain_str__name__[] = "constrain_str";

PyObject* pypdffit2_constrain_str(PyObject*, PyObject* args)
{
    PyObject* py_v = 0;
    char* vname;
    char* form;
    PyObject* py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OOss", &py_ppdf, &py_v, &vname, &form);
    if (!ok) return 0;
    PdfFit* ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    RefVar* v = (RefVar *) PyCapsule_GetPointer(py_v, cnvar);
    if (v->type() != "RefVar")
    {
        string emsg = "cannot constrain non-refinable variable ";
        emsg += vname;
        PyErr_SetString(pypdffit2_constraintError, emsg.c_str());
        return 0;
    }
    if (v->isAssigned()) {
        try {
            ppdf->constrain(*v, form);
        }
        catch (constraintError e) {
            PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
            return 0;
        }
        catch (unassignedError e) {
            PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
            return 0;
        }
    }
    else {
        ostringstream emsg;
        emsg << "Variable " << vname << " was not yet assigned";
        PyErr_SetString(pypdffit2_unassignedError, emsg.str().c_str());
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// constrain to integer.
char pypdffit2_constrain_int__doc__[] = "Constrain refined variable to integer.";
char pypdffit2_constrain_int__name__[] = "constrain_int";

PyObject* pypdffit2_constrain_int(PyObject*, PyObject* args)
{
    PyObject* py_v = 0;
    int ftype = 0;
    char* vname;
    int ipar;
    PyObject* py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OOsi|i", &py_ppdf, &py_v, &vname, &ipar, &ftype);
    if (!ok) return 0;
    PdfFit* ppdf = (PdfFit*) PyCapsule_GetPointer(py_ppdf, cnpfit);
    RefVar* v = (RefVar*) PyCapsule_GetPointer(py_v, cnvar);
    if (v->type() != "RefVar")
    {
        string emsg = "cannot constrain non-refinable variable ";
        emsg += vname;
        PyErr_SetString(pypdffit2_constraintError, emsg.c_str());
        return 0;
    }
    if (v->isAssigned()) {
        try {
            if (ftype)
            {
                ppdf->constrain(*v, ipar, (FCON)ftype);
            }
            else ppdf->constrain(*v, ipar);
        }
        catch (constraintError e) {
            PyErr_SetString(pypdffit2_constraintError, e.GetMsg().c_str());
            return 0;
        }
        catch (unassignedError e) {
            PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
            return 0;
        }
    }
    else {
        ostringstream emsg;
        emsg << "Variable " << vname << " was not yet assigned";
        PyErr_SetString(pypdffit2_unassignedError, emsg.str().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// setpar with value of double
char pypdffit2_setpar_dbl__doc__[] = "Set parameter value.";
char pypdffit2_setpar_dbl__name__[] = "setpar_dbl";

PyObject * pypdffit2_setpar_dbl(PyObject *, PyObject *args)
{
    unsigned int n;
    double val;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OId", &py_ppdf, &n, &val);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->setpar(n, val);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// setpar with value of RefVar
char pypdffit2_setpar_RV__doc__[] = "Set parameter value via refined variable.";
char pypdffit2_setpar_RV__name__[] = "setpar_RV";

PyObject * pypdffit2_setpar_RV(PyObject *, PyObject *args)
{
    unsigned int n;
    RefVar *v;
    PyObject *py_v;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OIO", &py_ppdf, &n, &py_v);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    v = (RefVar *) PyCapsule_GetPointer(py_v, cnvar);
    if( v->isAssigned() ) {
        try
        {
            ppdf->setpar(n, *v);
        }
        catch(unassignedError e)
        {
            PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
            return 0;
        }
    }
    else {
        string eout = "Variable not yet assigned";
        PyErr_SetString(pypdffit2_unassignedError, eout.c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// setvar
char pypdffit2_setvar__doc__[] = "Set variable to value.";
char pypdffit2_setvar__name__[] = "setvar";

PyObject * pypdffit2_setvar(PyObject *, PyObject *args)
{
    double a;
    NonRefVar *v;
    PyObject *py_v;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OOd", &py_ppdf, &py_v, &a);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    v = (NonRefVar *) PyCapsule_GetPointer(py_v, cnvar);
    if( v->isAssigned() ) {
        ppdf->setvar(*v, a);
    }
    else {
        string eout = "Must import a structure";
        PyErr_SetString(pypdffit2_unassignedError, eout.c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// getvar
char pypdffit2_getvar__doc__[] = "Get variable value.";
char pypdffit2_getvar__name__[] = "getvar";

PyObject * pypdffit2_getvar(PyObject *, PyObject *args)
{
    NonRefVar *v = 0;
    PyObject *py_v = 0;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OO", &py_ppdf, &py_v);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    v = (NonRefVar *) PyCapsule_GetPointer(py_v, cnvar);
    if(v->isAssigned()) {
        double crval = ppdf->getvar(*v);
        return Py_BuildValue("d", crval);
    }
    else {
        string eout = "Variable not yet assigned";
        PyErr_SetString(pypdffit2_unassignedError, eout.c_str());
        return 0;
        Py_INCREF(Py_None);
        return Py_None;
    }
}

// getcrw
char pypdffit2_getcrw__doc__[] = "Get cumulative Rw for the current dataset.";
char pypdffit2_getcrw__name__[] = "getcrw";

PyObject * pypdffit2_getcrw(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        vector<double> crw = ppdf->getcrw();
        PyObject *py_r;
        py_r = PyList_New(crw.size());
        for (int i = 0; i != int(crw.size()); ++i)
        {
            PyList_SetItem(py_r, i, PyFloat_FromDouble(crw[i]));
        }
        return py_r;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}

// getrw
char pypdffit2_getrw__doc__[] = "Get rw of fit.";
char pypdffit2_getrw__name__[] = "getrw";

PyObject * pypdffit2_getrw(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    double crval = ppdf->getrw();
    return Py_BuildValue("d", crval);
}

// getR
char pypdffit2_getR__doc__[] = "Get list of r-values for plotting.";
char pypdffit2_getR__name__[] = "getR";

PyObject * pypdffit2_getR(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    /* One should not put functionality in the bindings. However,
     * this function is meant to create a python object from a
     * c-object that does not actually exist. All that is stored
     * in pdffit about the fitted R-range is the minimum, maximum,
     * and step size. This binding turns that info into a python
     * list. If the number of fit points stored in nfmax and nfmin
     * have not yet been assigned, then this falls back on the
     * r-points. This is crucial so that one can investigate the
     * PDF before any fitting has been done. Once fit parameters
     * have been established, they can be used and the r list,
     * pdf-obs list, and pdf-fit list will be the same size.
     * Also see getpdf_obs().
     */
    try
    {
        int nfmin = ppdf->getnfmin();
        int nfmax = ppdf->getnfmax();
        int len = nfmax - nfmin + 1;
        double rmin = ppdf->getrmin();
        double rmax = ppdf->getrmax();
        double deltar = ppdf->getdeltar();
        PyObject *py_r;
        if(len == 1)
        {
            nfmin = 0;
            nfmax = (int) ((rmax - rmin)/deltar);
            len = nfmax + 1;
        }
        py_r = PyList_New(len);
        for (int i=nfmin;i<=nfmax;i++)
        {
            PyList_SetItem(py_r, i-nfmin, Py_BuildValue("d", i*deltar + rmin));
        }

        return py_r;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}


// getpdf_fit
char pypdffit2_getpdf_fit__doc__[] = "Get list of calculated pdf points.";
char pypdffit2_getpdf_fit__name__[] = "getpdf_fit";

PyObject * pypdffit2_getpdf_fit(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    //Return only the data range used in the fit
    try
    {
        int min = ppdf->getnfmin();
        int max = ppdf->getnfmax();
        int len = max - min + 1;
        vector<double> v_pdfdata = ppdf->getpdf_fit();
        PyObject *py_r;
        py_r = PyList_New(len);
        for (int i=min;i<=max;i++) {
            PyList_SetItem(py_r, i-min, Py_BuildValue("d", v_pdfdata[i]));
        }
        return py_r;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}


// getpdf_obs
char pypdffit2_getpdf_obs__doc__[] = "Get list of observed (theory) pdf points.";
char pypdffit2_getpdf_obs__name__[] = "getpdf_obs";

PyObject * pypdffit2_getpdf_obs(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        vector<double> v_pdfdata = ppdf->getpdf_obs();
        int nfmin = ppdf->getnfmin();
        int nfmax = ppdf->getnfmax();
        int len = nfmax - nfmin + 1;
        //Return only the data range used in the fit
        PyObject *py_r;
        py_r = PyList_New(len);
        for (int i=nfmin;i<=nfmax;i++) {
            PyList_SetItem(py_r, i-nfmin, Py_BuildValue("d", v_pdfdata[i]));
        }
        return py_r;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}

// getpdf_diff
char pypdffit2_getpdf_diff__doc__[] = "Get list of differences between observed and fitted PDF points.";
char pypdffit2_getpdf_diff__name__[] = "getpdf_diff";

PyObject * pypdffit2_getpdf_diff(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        vector<double> Gobs = ppdf->getpdf_obs();
        vector<double> Gfit = ppdf->getpdf_fit();
        int nfmin = ppdf->getnfmin();
        int nfmax = ppdf->getnfmax();
        int len = nfmax - nfmin + 1;
        //Return only the data range used in the fit
        PyObject *py_r;
        py_r = PyList_New(len);
        for (int i = nfmin; i <= nfmax; i++) {
            double Gdiff_i = Gobs[i] - Gfit[i];
            PyList_SetItem(py_r, i-nfmin, Py_BuildValue("d", Gdiff_i));
        }
        return py_r;
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}


// getpar
char pypdffit2_getpar__doc__[] = "Get value of parameter";
char pypdffit2_getpar__name__[] = "getpar";

PyObject * pypdffit2_getpar(PyObject *, PyObject *args)
{
    unsigned int n;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "OI", &py_ppdf, &n);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        double crval = ppdf->getpar(n);
        return Py_BuildValue("d", crval);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}

// fixpar
char pypdffit2_fixpar__doc__[] = "Fix value of parameter.";
char pypdffit2_fixpar__name__[] = "fixpar";

PyObject * pypdffit2_fixpar(PyObject *, PyObject *args)
{
    int n;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &n);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->fixpar(n);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// freepar
char pypdffit2_freepar__doc__[] = "Free parameter.";
char pypdffit2_freepar__name__[] = "freepar";

PyObject * pypdffit2_freepar(PyObject *, PyObject *args)
{
    int n;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &n);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->freepar(n);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// setphase
char pypdffit2_setphase__doc__[] = "Set phase in focus.";
char pypdffit2_setphase__name__[] = "setphase";

PyObject * pypdffit2_setphase(PyObject *, PyObject *args)
{
    int ip;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->setphase(ip);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// setdata
char pypdffit2_setdata__doc__[] = "Set data in focus.";
char pypdffit2_setdata__name__[] = "setdata";

PyObject * pypdffit2_setdata(PyObject *, PyObject *args)
{
    int is;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &is);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->setdata(is);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// psel
char pypdffit2_psel__doc__[] = "Select phase in focus.";
char pypdffit2_psel__name__[] = "psel";

PyObject * pypdffit2_psel(PyObject *, PyObject *args)
{
    int ip;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->selphase(ip);
    }
    catch(unassignedError e)
    {
       // PyErr_Warn(PyExc_Warning, e.GetMsg().c_str());
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// pdesel
char pypdffit2_pdesel__doc__[] = "Deselect phase in focus.";
char pypdffit2_pdesel__name__[] = "pdesel";

PyObject * pypdffit2_pdesel(PyObject *, PyObject *args)
{
    int ip;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->pdesel(ip);
    }
    catch(unassignedError e)
    {
       // PyErr_Warn(PyExc_Warning, e.GetMsg().c_str());
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// selectAtomType
char pypdffit2_selectAtomType__doc__[] = "Include element in 'i' or 'j' pair of PDF calculation.";
char pypdffit2_selectAtomType__name__[] = "selectAtomType";

PyObject * pypdffit2_selectAtomType(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip;
    char ijchar;
    char* smbpat;
    bool select;
    int ok = PyArg_ParseTuple(args, "Oicsb", &py_ppdf, &ip, &ijchar, &smbpat, &select);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->selectAtomType(ip, ijchar, smbpat, select);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e)
    {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// selectAtomIndex
char pypdffit2_selectAtomIndex__doc__[] = "Include atom of given index in 'i' or 'j' pair of PDF calculation.";
char pypdffit2_selectAtomIndex__name__[] = "selectAtomIndex";

PyObject * pypdffit2_selectAtomIndex(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip;
    char ijchar;
    int aidx1;
    bool select;
    int ok = PyArg_ParseTuple(args, "Oicib", &py_ppdf, &ip, &ijchar, &aidx1, &select);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->selectAtomIndex(ip, ijchar, aidx1, select);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e)
    {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// selectAll
char pypdffit2_selectAll__doc__[] = "Include all atoms in 'i' or 'j' pair of PDF calculation.";
char pypdffit2_selectAll__name__[] = "selectAll";

PyObject * pypdffit2_selectAll(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip;
    char ijchar;
    int ok = PyArg_ParseTuple(args, "Oic", &py_ppdf, &ip, &ijchar);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->selectAll(ip, ijchar);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e)
    {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// selectNone
char pypdffit2_selectNone__doc__[] = "Include all atoms in 'i' or 'j' pair of PDF calculation.";
char pypdffit2_selectNone__name__[] = "selectNone";

PyObject * pypdffit2_selectNone(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip;
    char ijchar;
    int ok = PyArg_ParseTuple(args, "Oic", &py_ppdf, &ip, &ijchar);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try
    {
        ppdf->selectNone(ip, ijchar);
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e)
    {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// bond_angle
char pypdffit2_bond_angle__doc__[] = "Return bond angle between three atoms.";
char pypdffit2_bond_angle__name__[] = "bond_angle";

PyObject * pypdffit2_bond_angle(PyObject *, PyObject *args)
{
    int ia, ja, ka;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oiii", &py_ppdf, &ia, &ja, &ka);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        pair<double,double> angstd = ppdf->bond_angle(ia, ja, ka);
        PyObject* py_tpl;
        py_tpl = Py_BuildValue("(d,d)", angstd.first, angstd.second);
        return py_tpl;
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    catch (unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
}

// bond_length_atoms (nearest bond length between two atoms)
char pypdffit2_bond_length_atoms__doc__[] =
    "Return a tuple of (dij, ddij) for distance between two atoms\n"
    "and its standard deviation.";
char pypdffit2_bond_length_atoms__name__[] = "bond_length_atoms";

PyObject * pypdffit2_bond_length_atoms(PyObject *, PyObject *args)
{
    int ia, ja;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oii", &py_ppdf, &ia, &ja);
    PairDistance pd;
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
	pd = ppdf->bond_length_atoms(ia, ja);
        PyObject *py_tpl;
        py_tpl = Py_BuildValue("(d,d)", pd.dij, pd.ddij);
	return py_tpl;
    }
    catch (ValueError e) {
	PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
	return 0;
    }
    catch (unassignedError e) {
	PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
	return 0;
    }
}

// bond_length_types (bond lengths between two elements inside given bounds)
char pypdffit2_bond_length_types__doc__[] =
    "Return bond lengths between two elements within given bounds\n"
    "\n"
    "a1    -- symbol of the first element in pair or 'ALL'\n"
    "a2    -- symbol of the second element in pair or 'ALL'\n"
    "lb    -- lower bound for bond lengths\n"
    "ub    -- upper bound for bond lengths\n"
    "\n"
    "Return a dictionary of distance data containing:\n"
    "\n"
    "dij  : list of bond lenghts within given bounds\n"
    "ddij : list of bond legnth standard deviations\n"
    "ij0  : list of tupled pairs of indices starting at 0\n"
    "ij1  : list of tupled pairs of indices starting at 1";
char pypdffit2_bond_length_types__name__[] = "bond_length_types";

PyObject * pypdffit2_bond_length_types(PyObject *, PyObject *args)
{
    char* symi;
    char* symj;
    double bmin = 0;
    double bmax = 0;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ossdd", &py_ppdf, &symi, &symj, &bmin, &bmax);
    vector<PairDistance> pdvec;
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
	pdvec = ppdf->bond_length_types(symi, symj, bmin, bmax);
        int np = pdvec.size();
        PyObject* py_dij;
        PyObject* py_ddij;
        PyObject* py_ij0;
        PyObject* py_ij1;
        py_dij = PyList_New(np);
        py_ddij = PyList_New(np);
        py_ij0 = PyList_New(np);
        py_ij1 = PyList_New(np);
	for (int i = 0; i < np; ++i)
	{
	    PairDistance& pd = pdvec[i];
	    PyObject *py_item;
            py_item = PyFloat_FromDouble(pd.dij);
            PyList_SetItem(py_dij, i, py_item);
            py_item = PyFloat_FromDouble(pd.ddij);
            PyList_SetItem(py_ddij, i, py_item);
	    py_item = Py_BuildValue("(i,i)", pd.i - 1, pd.j - 1);
            PyList_SetItem(py_ij0, i, py_item);
	    py_item = Py_BuildValue("(i,i)", pd.i, pd.j);
            PyList_SetItem(py_ij1, i, py_item);
        }
        PyObject* py_rv;
        py_rv = PyDict_New();
        PyDict_SetItemString(py_rv, "dij", py_dij);
        PyDict_SetItemString(py_rv, "ddij", py_ddij);
        PyDict_SetItemString(py_rv, "ij0", py_ij0);
        PyDict_SetItemString(py_rv, "ij1", py_ij1);
	return py_rv;
    }
    catch (ValueError e) {
	PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
	return 0;
    }
    catch (unassignedError e) {
	PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
	return 0;
    }
}

// get_scat_string
char pypdffit2_get_scat_string__doc__[] = "Return string with scatter details.";
char pypdffit2_get_scat_string__name__[] = "get_scat_string";

PyObject * pypdffit2_get_scat_string(PyObject *, PyObject *args)
{
    char stype;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oc", &py_ppdf, &stype);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    string outstring;
    if (!ppdf->curphase)
    {
        PyErr_SetString(pypdffit2_unassignedError, "No scatterers exist");
        return 0;
    }
    // here curphase exists, get_scat_string throws runtime error
    // for invalid stype.
    try {
        outstring = (ppdf->curphase)->get_scat_string(stype);
    }
    catch (runtime_error e) {
	PyErr_SetString(PyExc_ValueError, e.what());
        return 0;
    }
    return Py_BuildValue("s", outstring.c_str());
}

// get_scat
char pypdffit2_get_scat__doc__[] = "Return scattering factor for given element.";
char pypdffit2_get_scat__name__[] = "get_scat";

PyObject * pypdffit2_get_scat(PyObject *, PyObject *args)
{
    char stype;
    char* smbpat;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ocs", &py_ppdf, &stype, &smbpat);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    double value;
    try {
        value = ppdf->get_scat(stype, smbpat);
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    PyObject* py_rv;
    py_rv = PyFloat_FromDouble(value);
    return py_rv;
}

// set_scat
char pypdffit2_set_scat__doc__[] = "Set custom scattering factor for given element.";
char pypdffit2_set_scat__name__[] = "set_scat";

PyObject * pypdffit2_set_scat(PyObject *, PyObject *args)
{
    char stype;
    char* smbpat;
    double value;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Ocsd", &py_ppdf, &stype, &smbpat, &value);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    if (!ppdf->curphase)
    {
        PyErr_SetString(pypdffit2_unassignedError, "phase does not exist");
        return 0;
    }
    // Here curphase exists.  set_scat may throw
    // runtime or ValueError for invalid arguments
    try {
        ppdf->curphase->set_scat(stype, smbpat, value);
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    catch (runtime_error e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// reset_scat
char pypdffit2_reset_scat__doc__[] = "Reset scattering factor for given element.";
char pypdffit2_reset_scat__name__[] = "reset_scat";

PyObject * pypdffit2_reset_scat(PyObject *, PyObject *args)
{
    char* smbpat;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Os", &py_ppdf, &smbpat);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    if (!ppdf->curphase)
    {
        PyErr_SetString(pypdffit2_unassignedError, "phase does not exist");
        return 0;
    }
    // Here curphase exists.  reset_scat may throw
    // runtime or ValueError for invalid arguments
    try {
        ppdf->curphase->reset_scat(smbpat);
    }
    catch (ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

//throw-away variable that helps with exception handling
//RefVar *_junkRV = new RefVar();
// utility function (not available from Python)
RefVar *getRefVar(vector<RefVar> &v, unsigned int i)
{
    if (v.size() == 0)
    {
        throw unassignedError("Variable not yet assigned");
    }
    else if (i > v.size() || i < 1)
    {
        stringstream errstr;
        errstr << "Index " << i << " out of bounds";
        throw ValueError(errstr.str());
    }
    return &v[i-1];
}

//current phase and set refinable variable pointers

// lat
char pypdffit2_lat__doc__[] = "Pointer to refinable variable lat.";
char pypdffit2_lat__name__[] = "lat";

PyObject * pypdffit2_lat(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->lat,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// x
char pypdffit2_x__doc__[] = "Pointer to refinable variable x.";
char pypdffit2_x__name__[] = "x";

PyObject * pypdffit2_x(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->x,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// y
char pypdffit2_y__doc__[] = "Pointer to refinable variable y.";
char pypdffit2_y__name__[] = "y";

PyObject * pypdffit2_y(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->y,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// z
char pypdffit2_z__doc__[] = "Pointer to refinable variable z.";
char pypdffit2_z__name__[] = "z";

PyObject * pypdffit2_z(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->z,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u11
char pypdffit2_u11__doc__[] = "Pointer to refinable variable u11.";
char pypdffit2_u11__name__[] = "u11";

PyObject * pypdffit2_u11(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u11,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u22
char pypdffit2_u22__doc__[] = "Pointer to refinable variable u22.";
char pypdffit2_u22__name__[] = "u22";

PyObject * pypdffit2_u22(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u22,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u33
char pypdffit2_u33__doc__[] = "Pointer to refinable variable u33.";
char pypdffit2_u33__name__[] = "u33";

PyObject * pypdffit2_u33(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u33,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u12
char pypdffit2_u12__doc__[] = "Pointer to refinable variable u12.";
char pypdffit2_u12__name__[] = "u12";

PyObject * pypdffit2_u12(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u12,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u13
char pypdffit2_u13__doc__[] = "Pointer to refinable variable u13.";
char pypdffit2_u13__name__[] = "u13";

PyObject * pypdffit2_u13(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u13,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// u23
char pypdffit2_u23__doc__[] = "Pointer to refinable variable u23.";
char pypdffit2_u23__name__[] = "u23";

PyObject * pypdffit2_u23(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->u23,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// occ
char pypdffit2_occ__doc__[] = "Pointer to refinable variable occ.";
char pypdffit2_occ__name__[] = "occ";

PyObject * pypdffit2_occ(PyObject *, PyObject *args)
{
    int i;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "Oi", &py_ppdf, &i);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    try {
        RefVar *v = getRefVar(ppdf->occ,i);
        PyObject *py_v = PyCapsule_New(v, cnvar, NULL);
        return py_v;
    }
    catch(unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    catch(ValueError e) {
        PyErr_SetString(PyExc_ValueError, e.GetMsg().c_str());
        return 0;
    }
}

// pscale
char pypdffit2_pscale__doc__[] = "Pointer to variable pscale.";
char pypdffit2_pscale__name__[] = "pscale";

PyObject * pypdffit2_pscale(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->pscale), cnvar, NULL);
    return py_v;
}

// spdiameter
char pypdffit2_spdiameter__doc__[] = "Pointer to variable spdiameter.";
char pypdffit2_spdiameter__name__[] = "spdiameter";

PyObject * pypdffit2_spdiameter(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->spdiameter), cnvar, NULL);
    return py_v;
}

// stepcut
char pypdffit2_stepcut__doc__[] = "Pointer to nonvariable stepcut.";
char pypdffit2_stepcut__name__[] = "stepcut";

PyObject * pypdffit2_stepcut(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->stepcut), cnvar, NULL);
    return py_v;
}

// sratio
char pypdffit2_sratio__doc__[] = "Pointer to variable sratio.";
char pypdffit2_sratio__name__[] = "sratio";

PyObject * pypdffit2_sratio(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->sratio), cnvar, NULL);
    return py_v;
}

// delta2
char pypdffit2_delta2__doc__[] = "Pointer to variable delta2.";
char pypdffit2_delta2__name__[] = "delta2";

PyObject * pypdffit2_delta2(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->delta2), cnvar, NULL);
    return py_v;
}

// delta1
char pypdffit2_delta1__doc__[] = "Pointer to variable delta1.";
char pypdffit2_delta1__name__[] = "delta1";

PyObject * pypdffit2_delta1(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->delta1), cnvar, NULL);
    return py_v;
}

// dscale
char pypdffit2_dscale__doc__[] = "Pointer to variable dscale.";
char pypdffit2_dscale__name__[] = "dscale";

PyObject * pypdffit2_dscale(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->dscale), cnvar, NULL);
    return py_v;
}

// qdamp
char pypdffit2_qdamp__doc__[] = "Pointer to variable qdamp.";
char pypdffit2_qdamp__name__[] = "qdamp";

PyObject * pypdffit2_qdamp(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->qdamp), cnvar, NULL);
    return py_v;
}

// qbroad
char pypdffit2_qbroad__doc__[] = "Pointer to variable qbroad.";
char pypdffit2_qbroad__name__[] = "qbroad";

PyObject * pypdffit2_qbroad(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->qbroad), cnvar, NULL);
    return py_v;
}

// rcut
char pypdffit2_rcut__doc__[] = "Pointer to nonvariable rcut.";
char pypdffit2_rcut__name__[] = "rcut";

PyObject * pypdffit2_rcut(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    PyObject *py_v = PyCapsule_New(&(ppdf->rcut), cnvar, NULL);
    return py_v;
}

// get_atoms
char pypdffit2_get_atoms__doc__[] = "Get element symbols of atoms in the phase.";
char pypdffit2_get_atoms__name__[] = "get_atoms";

PyObject * pypdffit2_get_atoms(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip = 0;
    int ok = PyArg_ParseTuple(args, "O|i", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    Phase* ph;
    try {
	ph = ppdf->getphase(ip);
    }
    catch (unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    // Phase ph is defined here
    PyObject *py_atoms = PyList_New(ph->natoms);
    for (int i = 0; i < ph->natoms; ++i)
    {
	string usymbol = toupper(ph->atom[i].atom_type->symbol);
        PyList_SetItem(py_atoms, i, PyUnicode_FromString(usymbol.c_str()));
    }
    return py_atoms;
}

// num_atoms
char pypdffit2_num_atoms__doc__[] = "Get the number of atoms in the current phase.";
char pypdffit2_num_atoms__name__[] = "num_atoms";

PyObject * pypdffit2_num_atoms(PyObject *, PyObject *args)
{
    int retval = 0;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    if (ppdf->curphase)
    {
        retval = (ppdf->curphase)->natoms;
    }
    else
    {
        PyErr_SetString(pypdffit2_unassignedError, "No data loaded");
        return 0;
    }
    return Py_BuildValue("i", retval);
}

// get_atom_types
char pypdffit2_get_atom_types__doc__[] = "Get ordered unique symbols of atoms in the phase.";
char pypdffit2_get_atom_types__name__[] = "get_atom_types";

PyObject * pypdffit2_get_atom_types(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ip = 0;
    int ok = PyArg_ParseTuple(args, "O|i", &py_ppdf, &ip);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    Phase* ph;
    try {
	ph = ppdf->getphase(ip);
    }
    catch (unassignedError e) {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    // Phase ph is defined here
    PyObject *py_atom_types = PyList_New(ph->atom_types.size());
    for (int i = 0; i < int(ph->atom_types.size()); ++i)
    {
	string usymbol = toupper(ph->atom_types[i]->symbol);
        PyList_SetItem(py_atom_types, i, PyUnicode_FromString(usymbol.c_str()));
    }
    return py_atom_types;
}


// num_phases
char pypdffit2_num_phases__doc__[] = "Get the number of loaded phases.";
char pypdffit2_num_phases__name__[] = "num_phases";

PyObject * pypdffit2_num_phases(PyObject *, PyObject *args)
{
    int retval = 0;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    retval = ppdf->num_phases();
    return Py_BuildValue("i", retval);
}


// num_datasets
char pypdffit2_num_datasets__doc__[] = "Get the number of loaded datasets.";
char pypdffit2_num_datasets__name__[] = "num_datasets";

PyObject * pypdffit2_num_datasets(PyObject *, PyObject *args)
{
    int retval = 0;
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    retval = ppdf->num_datasets();
    return Py_BuildValue("i", retval);
}


// phase_fractions
char pypdffit2_phase_fractions__doc__[] =
    "Return relative phase fractions for current dataset scattering type\n"
    "\n"
    "Return a dictionary of relative phase fractions:\n"
    "\n"
    "atom    -- list of fractions normalized to atom count\n"
    "stdatom -- errors of atom count fractions\n"
    "cell    -- list of fractions normalized to unit cell count\n"
    "stdcell -- errors of unit cell count fractions\n"
    "mass    -- list of relative weight fractions\n"
    "stdmass -- errors of relative weight fractions\n"
    ;
char pypdffit2_phase_fractions__name__[] = "phase_fractions";

PyObject * pypdffit2_phase_fractions(PyObject *, PyObject *args)
{
    PyObject *py_ppdf = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_ppdf);
    if (!ok) return 0;
    PdfFit *ppdf = (PdfFit *) PyCapsule_GetPointer(py_ppdf, cnpfit);
    map <string, vector<double> > fractions;
    try
    {
        fractions = ppdf->getPhaseFractions();
    }
    catch(unassignedError e)
    {
        PyErr_SetString(pypdffit2_unassignedError, e.GetMsg().c_str());
        return 0;
    }
    // convert fractions map to an equivalent python dictionary
    map <string, vector<double> >::iterator ii;
    PyObject* py_rv;
    py_rv = PyDict_New();
    for (ii = fractions.begin(); ii != fractions.end(); ++ii)
    {
        int n = ii->second.size();
        PyObject* py_lst;
        py_lst = PyList_New(n);
        for (int i = 0; i < n; ++i)
        {
            PyObject* py_value;
            py_value = PyFloat_FromDouble(ii->second.at(i));
            PyList_SetItem(py_lst, i, py_value);
        }
        PyDict_SetItemString(py_rv, ii->first.c_str(), py_lst);
    }
    return py_rv;
}

// redirect_stdout
char pypdffit2_redirect_stdout__doc__[] = "Redirect engine output to a file-like object.";
char pypdffit2_redirect_stdout__name__[] = "redirect_stdout";

PyObject * pypdffit2_redirect_stdout(PyObject *, PyObject *args)
{
    // instance of PyFileStreambuf which takes care of redirection
    PyObject *py_file = 0;
    int ok = PyArg_ParseTuple(args, "O", &py_file);
    if (!ok) return 0;
    // check if py_file has write and flush attributes
    if ( !PyObject_HasAttrString(py_file, "write") ||
	 !PyObject_HasAttrString(py_file, "flush") )
    {
        PyErr_SetString(PyExc_TypeError, "expected file-like argument");
        return 0;
    }
    // create py_stdout_streambuf if necessary
    if (!py_stdout_streambuf)
    {
	py_stdout_streambuf = new PyFileStreambuf(py_file);
    }
    py_stdout_streambuf->redirect(py_file);
    // on first redirection we need to assign new ostream to NS_PDFFIT2::pout
    if (NS_PDFFIT2::pout == &std::cout)
    {
	NS_PDFFIT2::pout = new ostream(py_stdout_streambuf);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// restore_stdout
char pypdffit2_restore_stdout__doc__[] =
    "Restore engine output to the default stream (std::cout).";
char pypdffit2_restore_stdout__name__[] =
    "restore_stdout";

PyObject * pypdffit2_restore_stdout(PyObject *, PyObject *args)
{
    // no arguments.
    if (!PyArg_ParseTuple(args, ""))
        return 0;

    // If the global output stream pointer is not std::cout, then delete the custom stream.
    if (NS_PDFFIT2::pout != &std::cout)
    {
        delete NS_PDFFIT2::pout;
        NS_PDFFIT2::pout = &std::cout;
    }

    // Clean up the custom stream buffer
    if (py_stdout_streambuf)
    {
        delete py_stdout_streambuf;
        py_stdout_streambuf = nullptr;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

// is_element
char pypdffit2_is_element__doc__[] = "Check if element or isotope is defined in the built-in periodic table.";
char pypdffit2_is_element__name__[] = "is_element";

PyObject * pypdffit2_is_element(PyObject *, PyObject *args)
{
    // instance of PyFileStreambuf which takes care of redirection
    char *smbpat;
    int ok = PyArg_ParseTuple(args, "s", &smbpat);
    if (!ok) return 0;
    const LocalPeriodicTable* lpt = LocalPeriodicTable::instance();
    PyObject *rv = PyBool_FromLong(lpt->has(smbpat));
    return rv;
}

// End of file
