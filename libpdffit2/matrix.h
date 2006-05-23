#if !defined(INCL_MATRIX)

#define INCL_MATRIX

#include <vector>
using namespace std;

template <class Type> class matrix
{
	vector<vector<Type> > a;
	int rows, cols;

	public:
	matrix() {}
	matrix(int m,int n) { resize(m,n); }
	void resize(int m, int n) { a.resize(m); for (int i=0; i<m; i++) a[i].resize(n);  rows=m; cols=n;}
	void clear() { for (unsigned int i=0; i<a.size(); i++) a[i].clear(); a.clear(); rows = cols = 0; }
	vector<Type>& operator[](int i) { return a[i]; }
	int getrows() { return rows; }
	int getcols() { return cols; }
	template <class T> friend ostream &operator<<(ostream &ostream, matrix<T> &v);

	matrix<Type> transposed() 
	{ 
		matrix<Type> &a=*this, b(a.cols,a.rows);

		for (int i=0; i<a.rows; i++)
			for (int j=0; j<a.cols; j++)
				b[j][i] = a[i][j];
		return b;
	}
	
	vector<Type> sd()   // returns standard deviation vector from covariance matrix
	{
		matrix<Type> &a=*this;
		vector<Type> v(a.cols);
		
		if (a.cols != a.rows)
		{
			cout << "Matrix not square in <diagonal>\n";
			throw "Matrix not square in <diagonal>";
		}
		
		for (int i=0; i<a.rows; i++)
			v[i] = sqrt(a[i][i]);
		
		return v;
	}

	
    matrix<Type> operator*(matrix<Type> &b)
    {
		matrix<Type> &a=*this;
        matrix<Type> c(a.rows,b.cols);
		
		if (a.cols != b.rows)
		{
			cout << "Inconsistent matrix multiplication\n";
			throw "Inconsistent matrix multiplication";
		}

		for (int i=0; i<a.rows; i++)
		{
			for (int j=0; j<b.cols; j++)
			{
				c[i][j] = 0;
				
				for (int k=0; k< a.cols; k++)
					c[i][j] += a[i][k]*b[k][j];
			}
		}
        return c;
    }
	
	void print() 
	{
		for (int i=0; i<rows; i++)
		{
			cout << i << ": ";
			for (int j=0; j<cols; j++)
				cout << a[i][j] << " ";
			cout << endl;
		}
	}
};

template <class Type> class tensor
{
	vector<matrix<Type> > a;

	public:
	tensor() {}
	tensor(int m,int n,int q) { a.resize(m); for (int i=0; i<m; i++) a[i].matrix(n,q); } 
	matrix<Type>& operator[](int i) { return a[i]; }
};

template <class Type> ostream &operator<<(ostream &ostream, vector<Type> &v)
{
	ostream << "(";
	for (int i=0; i<v.size(); i++) 
	{
		ostream << v[i];
		if (i!=(v.size()-1)) ostream << ", ";
		else ostream << ")\n";
	}
	return ostream;
}

template <class Type> ostream &operator<<(ostream &ostream, matrix<Type> &v)
{
	ostream << "( ";
	for (int i=0; i<v.rows; i++)
	{
		ostream << "(";
		for (int j=0; j<v.cols; j++)
		{
			cout << v.a[i][j];
			if (j!=(v.cols-1)) cout << ", "; 
			else cout << ")";
		}
		if (i!=(v.rows-1))  cout << ", ";
		else cout << " )\n";
	}
	return ostream;
}

#endif


