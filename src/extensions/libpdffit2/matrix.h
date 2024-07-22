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
* template class matrix
*
* Comments: optimized from original vector of vectors
*
***********************************************************************/

#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <algorithm>
#include <vector>

#include "OutputStreams.h"
using NS_PDFFIT2::pout;

template <class T> class matrix
{
    private:

	// Data Methods
	T* mdata;
	size_t mrows;
	size_t mcols;
	size_t msize;

    public:

	// Constructors
	matrix() :
	    mdata(NULL), mrows(0), mcols(0), msize(0)
	{ }
	matrix(const matrix& src) :
	    mdata(new T[src.mrows*src.mcols]),
	    mrows(src.mrows), mcols(src.mcols), msize(src.msize)
	{
	    std::copy(src.mdata, src.mdata + src.msize, mdata);
	}
	matrix(size_t m, size_t n) :
	    mrows(m), mcols(n), msize(mrows*mcols)
	{
	    mdata = new T[mrows*mcols];
	    std::fill(mdata, mdata + msize, T(0));
	}
	// Destructor
	~matrix()
	{
	    delete[] mdata;
	}
	// Methods
	matrix& operator=(const matrix& src)
	{
	    if (this == &src)	return *this;
	    delete[] mdata;
	    mdata = new T[src.msize];
	    std::copy(src.mdata, src.mdata + src.msize, mdata);
	    mrows = src.mrows;
	    mcols = src.mcols;
	    msize = src.msize;
	    return *this;
	}
	matrix& operator=(T value)
	{
	    std::fill_n(mdata, msize, value);
	    return *this;
	}
	void resize(size_t m, size_t n, T value=T(0))
	{
	    if (m == mrows && n == mcols)   return;
	    T* resized = new T[m*n];
	    std::fill(resized, resized + m*n, value);
	    for (   size_t i = 0, offset = 0;
		    i != std::min(m, mrows); ++i, offset += n )
	    {
		for (size_t j = 0; j != std::min(n, mcols); ++j)
		{
		    resized[offset+j] = (*this)(i, j);
		}
	    }
	    delete[] mdata;
	    mdata = resized;
	    mrows = m;
	    mcols = n;
	    msize = m*n;
	}
	void clear()
	{
	    delete[] mdata;
	    mdata = NULL;
	    mrows = mcols = msize = 0;
	}
	std::vector<T> rowVector(size_t i)
	{
	    return std::vector<T>(mdata + mcols*i, mdata + mcols*(i+1));
	}
	std::vector<T> colVector(size_t j)
	{
	    std::vector<T> column(mrows);
	    typename std::vector<T>::iterator vii;
	    vii = column.begin();
	    T* pt = mdata + j;
	    for (; vii != column.end(); ++vii, pt += mcols)   *vii = *pt;
	    return column;
	}
	inline size_t getrows()
	{
	    return mrows;
	}
	inline size_t getcols()
	{
	    return mcols;
	}
	inline T* operator[](size_t i)
	{
	    return mdata + mcols*i;
	}
	inline T& operator()(size_t i, size_t j)
	{
	    return *(mdata + i*mcols + j);
	}
	matrix<T> transposed()
	{
	    matrix<T> mxt(mcols, mrows);
	    T* pt = mxt.mdata;
	    for (size_t j = 0; j != mcols; ++j)
	    {
		for (size_t i = 0; i != mrows; ++i, ++pt)
		{
		    *pt = (*this)(i, j);
		}
	    }
	    return mxt;
	}
	// returns standard deviation vector from covariance matrix
	std::vector<T> sd()
	{
	    std::vector<T> v(mcols);
	    if (mcols != mrows)
	    {
		const char* emsg = "Matrix not square in <diagonal>";
		throw emsg;
	    }
	    typename std::vector<T>::iterator vii;
	    vii = v.begin();
	    for (T* pii = mdata; pii < mdata + msize; pii += mcols + 1, ++vii)
	    {
		*vii = sqrt(*pii);
	    }
	    return v;
	}
	matrix<T> operator*(matrix<T>& B)
	{
	    matrix<T>& A = *this;
	    matrix<T> C(A.mrows, B.mcols);
	    if (A.mcols != B.mrows)
	    {
		const char* emsg = "Inconsistent matrix multiplication";
		throw emsg;
	    }
	    for (size_t i = 0; i != A.mrows; ++i)
	    {
		for (size_t j = 0; j != B.mcols; ++j)
		{
		    T& Cij = C(i,j);
		    Cij = 0;
		    for (size_t k = 0; k != A.mcols; ++k)
		    {
			Cij += A(i,k)*B(k,j);
		    }
		}
	    }
	    return C;
	}
	template <class T1>
	    friend std::ostream& operator<<(std::ostream &out, matrix<T1> &mx);
	void print()
	{
	    using namespace std;
	    for (size_t i = 0; i != mrows; ++i)
	    {
		*pout << i << ": ";
		for (size_t j = 0; j != mcols; ++j)
		{
		    *pout << (*this)(i,j) << " ";
		}
		*pout << endl;
	    }
	}
};

template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> &v)
{
    using namespace std;
    out << "(";
    for (size_t i=0; i != v.size(); ++i)
    {
	out << v[i];
	if (i != (v.size()-1))	out << ", ";
	else			out << ")\n";
    }
    return out;
}

template <class T>
std::ostream& operator<<(std::ostream &out, matrix<T> &mx)
{
    using namespace std;
    out << "( ";
    for (size_t i = 0; i != mx.mrows; ++i)
    {
	out << "(";
	for (size_t j = 0; j != mx.mcols; ++j)
	{
	    *pout << mx(i,j);
	    if (j != (mx.mcols-1))  *pout << ", ";
	    else		    *pout << ")";
	}
	if (i != (mx.mrows-1))	    *pout << ", ";
	else			    *pout << " )\n";
    }
    return out;
}

#endif	// MATRIX_H_INCLUDED
