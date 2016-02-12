/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Utilities from numerical recipies.
*
* Comments:
*
***********************************************************************/

#include <cstdio>
#include <cstdlib>

const int getNR_END()
{
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
    /* allocate a vector with subscript range v[nl..nh] */
{
    T *v = NULL;
    if (nl > nh)  return v;
    v=(T *)malloc((size_t) ((nh-nl+1+getNR_END())*sizeof(T)));
    if (!v) nrerror("allocation failure in _vector()");
    return v-nl+getNR_END();
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

template <class T> void _free_vector(T *v, long nl, long nh)
    /* free a <class T> vector allocated with vector() */
{
    if (nl > nh)  return;
    free((T*) (v+nl-getNR_END()));
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

template <class T> T **_matrix(long nrl, long nrh, long ncl, long nch)
    /* allocate a <class T> matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    T **m = NULL;
    if (nrl > nrh || ncl > nch)  return m;

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

double **dmatrix(long nrl, long nrh, long ncl, long nch)
    /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    return _matrix<double>(nrl, nrh, ncl, nch);
}

template <class T> void _free_matrix(T **m, long nrl, long nrh, long ncl, long nch)
    /* free a double matrix allocated by matrix() */
{
    if (nrl > nrh || ncl > nch)  return;
    free((T*) (m[nrl]+ncl-getNR_END()));
    free((T**) (m+nrl-getNR_END()));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
    /* free a double matrix allocated by matrix() */
{
    _free_matrix<double>(m, nrl, nrh, ncl, nch);
}

// End of file
