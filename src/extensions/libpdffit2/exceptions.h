/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Chris Farrow, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Exceptions used in pdffit2.
*
* Comments:
*
***********************************************************************/

#ifndef EXCEPTIONS_H_INCLUDED
#define EXCEPTIONS_H_INCLUDED

#include <iostream>
#include <string>

class Exception
{
    private:
	std::string msg;
    public:
	Exception(std::string _msg) : msg(_msg) {}
	std::string GetMsg()
	{
	    return msg;
	}
};

//specific exceptions - mimic python names
class ValueError : public Exception
{
    public:
	ValueError(std::string _msg) : Exception(_msg) {}
};

class unassignedError : public Exception
{
    public:
	unassignedError(std::string _msg) : Exception(_msg) {}
};

class IOError : public Exception
{
    public:
	IOError(std::string _msg) : Exception(_msg) {}
};

class dataError : public Exception
{
    public:
	dataError(std::string _msg) : Exception(_msg) {}
};

class structureError : public Exception
{
    public:
	structureError(std::string _msg) : Exception(_msg) {}
};

class constraintError : public Exception
{
    public:
	constraintError(std::string _msg) : Exception(_msg) {}
};

class calculationError : public Exception
{
    public:
	calculationError(std::string _msg) : Exception(_msg) {}
};

class parseError : public Exception
{
    public:
	parseError(std::string _msg) : Exception(_msg) {}
};

//This one is used internally, and should not make it to the python layer.
class vgetException : public Exception
{
    public:
	vgetException(std::string _msg) : Exception(_msg) {}
};

#endif	// EXCEPTIONS_H_INCLUDED
