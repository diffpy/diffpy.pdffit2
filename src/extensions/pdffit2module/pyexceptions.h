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

#ifndef PYPDFFIT2_EXCEPTIONS_H
#define PYPDFFIT2_EXCEPTIONS_H

// exceptions

extern PyObject *pypdffit2_runtimeError;
extern PyObject *pypdffit2_unassignedError;
extern PyObject *pypdffit2_dataError;
extern PyObject *pypdffit2_structureError;
extern PyObject *pypdffit2_constraintError;
extern PyObject *pypdffit2_calculationError;

#endif	// PYPDFFIT2_EXCEPTIONS_H
