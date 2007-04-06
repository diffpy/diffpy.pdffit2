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
* The python pdffit2 module.
*
* Comments:
*
* $Id$
*
***********************************************************************/

#include <Python.h>
#include <ostream>

#include "pyexceptions.h"
#include "bindings.h"
#include "libpdffit2/OutputStreams.h"

using namespace std;

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

// End of file
