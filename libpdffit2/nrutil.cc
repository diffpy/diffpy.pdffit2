#include <stdio.h>
#include <stdlib.h>

//#define NR_END 1
//#define FREE_ARG char*

const int getNR_END() {
    return 1;
}

static void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

template <class T> T *_vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	T *v;
	v=(T *)malloc((size_t) ((nh-nl+1+getNR_END())*sizeof(T)));
	if (!v) nrerror("allocation failure in _vector()");
	return v-nl+getNR_END();
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	return _vector<float>(nl, nh);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	return _vector<double>(nl, nh);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	return _vector<int>(nl, nh);
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	return _vector<unsigned long>(nl, nh);
}

template <class T> void _free_vector(T *v, long nl, long nh)
/* free a <class T> vector allocated with vector() */
{
	free((char*) (v+nl-getNR_END()));
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	_free_vector<float>(v, nl, nh);
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	_free_vector<double>(v, nl, nh);
}


void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	_free_vector<int>(v, nl, nh);
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	_free_vector<unsigned long>(v, nl, nh);
}

template <class T> T **_matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a <class T> matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	T **m;

	/* allocate pointers to rows */
	m=(T **) malloc((size_t)((nrow+getNR_END())*sizeof(T*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += getNR_END();
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(T *) malloc((size_t)((nrow*ncol+getNR_END())*sizeof(T)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += getNR_END();
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	return _matrix<float>(nrl, nrh, ncl, nch);
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	return _matrix<double>(nrl, nrh, ncl, nch);
}

template <class T> void _free_matrix(T **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((char*) (m[nrl]+ncl-getNR_END()));
	free((char*) (m+nrl-getNR_END()));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	_free_matrix<float>(m, nrl, nrh, ncl, nch);
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	_free_matrix<double>(m, nrl, nrh, ncl, nch);
}
