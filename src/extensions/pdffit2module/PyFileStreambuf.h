/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class PyFileStreambuf
*
* Comments: PyFileStreambuf is a C++ streambuf which writes to a python
*           file-like object.  The python file can be changed anytime by
*           calling the redirect() method.
*
* Examples: // redirect std::cout
*           std::cout.rdbuf( PyFileStreambuf(PyObject* python_file) );
*
***********************************************************************/

#ifndef PYFILESTREAMBUF_H_INCLUDED
#define PYFILESTREAMBUF_H_INCLUDED
#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <streambuf>

// MS compatibility fix
#include <memory>

class PyFileStreambuf : public std::streambuf
{
    private:

        // Data members
        PyObject* py_file;

    public:

        // Constructor
        PyFileStreambuf(PyObject* f) : py_file(f)
        {
            Py_INCREF(py_file);
        }

        // Destructor
        ~PyFileStreambuf()
        {
            Py_DECREF(py_file);
        }

        // Methods
        PyObject* redirect(PyObject* f)
        {
            Py_INCREF(f);
            Py_DECREF(py_file);
            py_file = f;
            return py_file;
        }

    protected:

        virtual int_type overflow( int_type c)
        {
            char ch = static_cast<char>(c);
            PyObject* py_str = PyUnicode_FromStringAndSize(&ch, 1);
            if (!py_str) { return traits_type::eof(); }
            PyObject* rv = PyObject_CallMethod(py_file, "write", "O", py_str);
            Py_DECREF(py_str);
            if (rv)  { Py_DECREF(rv); }
            return c;
        }

        virtual std::streamsize xsputn(const char_type* s, std::streamsize n)
        {
            PyObject* py_str = PyUnicode_DecodeUTF8(s, n, "replace");
            if (!py_str) { return 0; }
            PyObject* rv = PyObject_CallMethod(py_file, "write", "O", py_str);
            Py_DECREF(py_str);
            if (rv)  { Py_DECREF(rv); }
            return n;
        }

        virtual int sync()
        {
            if (PyObject_HasAttrString(py_file, "flush"))
            {
                PyObject* rv;
                rv = PyObject_CallMethod(py_file, "flush", NULL);
                if (rv)  { Py_DECREF(rv); }
            }
            return 0;
        }

};

#endif  // PYFILESTREAMBUF_H_INCLUDED
