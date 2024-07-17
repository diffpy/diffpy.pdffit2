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
* Exceptions for pdffit2 python module.
*
* Comments:
*
***********************************************************************/
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *pypdffit2_runtimeError = 0;
PyObject *pypdffit2_unassignedError = 0;
PyObject *pypdffit2_dataError = 0;
PyObject *pypdffit2_structureError = 0;
PyObject *pypdffit2_calculationError = 0;
PyObject *pypdffit2_constraintError = 0;

// End of file
