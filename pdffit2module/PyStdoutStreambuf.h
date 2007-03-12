/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class PyStdoutStreambuf
*
* Comments: PyStdoutStreambuf is a C++ streambuf which writes to python
*	    sys.stdout.  PyStdoutStreambuf is a singleton class, the
*	    pointer can be obtained via PyStdoutStreambuf::instance().
*	    Intended for redirecting output to sys.stdout.
*
* Examples: // redirect std::cout
*	    std::cout.rdbuf( PyStdoutStreambuf::instance() );
*
*	    // create new output stream linked to sys.stdout
*	    std::ostream sys_stdout( PyStdoutStreambuf::instance() );
*
* $Id$
*
***********************************************************************/

#ifndef PYSTDOUTSTREAMBUF_H_INCLUDED
#define PYSTDOUTSTREAMBUF_H_INCLUDED

#include <Python.h>
#include <streambuf>

// MS compatibility fix
#include <memory>

class PyStdoutStreambuf : public std::streambuf
{
    private:

	// Data members
	PyObject* sys;

	// PyStdoutStreambuf is a singleton class
	PyStdoutStreambuf() :
	    sys(PyImport_ImportModule("sys"))
	{ }

	// sys.stdout
	PyObject* py_file()
	{
	    PyObject* out = PyObject_GetAttrString(sys, "stdout");
	    return out;
	}

    public:

	// Access to singleton instance
	static PyStdoutStreambuf* instance()
	{
	    static std::auto_ptr<PyStdoutStreambuf> p(new PyStdoutStreambuf());
	    return p.get();
	}
	// Destructor
	~PyStdoutStreambuf()
	{
	    Py_DECREF(sys);
	}

    protected:

	virtual int_type overflow( int_type c)
	{
	    PyObject* rv;
	    PyObject* out = py_file();
	    rv = PyObject_CallMethod(out, "write", "(s#)", &c, 1);
	    if (rv)	    Py_DECREF(rv);
	    Py_DECREF(out);
	    return c;
	}

	virtual std::streamsize xsputn(const char_type* s, std::streamsize n)
	{
	    PyObject* rv;
	    PyObject* out = py_file();
	    rv = PyObject_CallMethod(out, "write", "(s#)", s, n);
	    if (rv)	    Py_DECREF(rv);
	    Py_DECREF(out);
	    return n;
	}

	virtual int sync()
	{
	    PyObject* out = py_file();
	    if (PyObject_HasAttrString(out, "flush"))
	    {
		PyObject* rv;
		rv = PyObject_CallMethod(out, "flush", NULL);
		if (rv)	    Py_DECREF(rv);
	    }
	    Py_DECREF(out);
	    return 0;
	}

};

#endif	// PYSTDOUTSTREAMBUF_H_INCLUDED
