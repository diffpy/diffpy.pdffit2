// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Michael A.G. Aivazis
//                        California Institute of Technology
//                        (C) 1998-2005  All Rights Reserved
// 
//  <LicenseText>
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>

#include <Python.h>
//#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_pdffit2
//#include "Numeric/arrayobject.h"

#include "exceptions.h"
#include "bindings.h"


char pypdffit2_module__doc__[] = "Pdffit2";

// Initialization function for the module (*must* be called initpdffit2)
extern "C"
void
initpdffit2()
{
    // create the module and add the functions
    PyObject * m = Py_InitModule4(
        "pdffit2", pypdffit2_methods,
        pypdffit2_module__doc__, 0, PYTHON_API_VERSION);

    // make Numeric functions available
    //import_array();

    // get its dictionary
    PyObject * d = PyModule_GetDict(m);

    // check for errors
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module pdffit2");
    }

    // install the module exceptions
    pypdffit2_runtimeError = PyErr_NewException("pdffit2.runtime", 0, 0);
    PyDict_SetItemString(d, "RuntimeException", pypdffit2_runtimeError);

    pypdffit2_unassignedError = PyErr_NewException(
            "pdffit2.unassignedError", 0, 0);
    PyDict_SetItemString(d, "unassignedError", pypdffit2_unassignedError);

    pypdffit2_dataError = PyErr_NewException(
            "pdffit2.dataError", 0, 0);
    PyDict_SetItemString(d, "dataError", pypdffit2_dataError);

    pypdffit2_structureError = PyErr_NewException(
            "pdffit2.structureError", 0, 0);
    PyDict_SetItemString(d, "structureError", pypdffit2_structureError);

    pypdffit2_calculationError = PyErr_NewException(
            "pdffit2.calculationError", 0, 0);
    PyDict_SetItemString(d, "calculationError", pypdffit2_calculationError);

    pypdffit2_constraintError = PyErr_NewException(
            "pdffit2.constraintError", 0, 0);
    PyDict_SetItemString(d, "constraintError", pypdffit2_constraintError);


    return;
}

// version
// $Id$

// End of file
