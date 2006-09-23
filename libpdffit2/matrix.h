#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <vector>

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
	void resize(size_t m, size_t n)
	{
	    if (m == mrows && n == mcols)   return;
	    T* resized = new T[m*n];
	    std::fill(resized, resized + m*n, T(0));
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
	std::vector<T> sd()   // returns standard deviation vector from covariance matrix
	{
	    using namespace std;
	    vector<T> v(mcols);
	    if (mcols != mrows)
	    {
		const char* emsg = "Matrix not square in <diagonal>";
		cerr << emsg << endl;
		throw emsg;
	    }
	    typename vector<T>::iterator vii;
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
		std::cerr << emsg << std::endl;
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
		cout << i << ": ";
		for (size_t j = 0; j != mcols; ++j)
		{
		    cout << (*this)(i,j) << " ";
		}
		cout << endl;
	    }
	}
};

/* not used anywhere
template <class Type> class tensor
{
	vector<matrix<Type> > a;

	public:
	tensor() {}
	tensor(int m,int n,int q) { a.resize(m); for (int i=0; i<m; i++) a[i].matrix(n,q); } 
	matrix<Type>& operator[](int i) { return a[i]; }
};
*/

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
	    cout << mx(i,j);
	    if (j != (mx.mcols-1))  cout << ", "; 
	    else		    cout << ")";
	}
	if (i != (mx.mrows-1))	    cout << ", ";
	else			    cout << " )\n";
    }
    return out;
}

#endif	// MATRIX_H_INCLUDED
