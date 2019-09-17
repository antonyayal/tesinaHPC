#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

float *matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i,j,N,
        nrow=nrh-nrl+1, ncol=nch-ncl+1;
    N = (nrow+1)*(ncol+1);
    float *m;

    m=(float *) malloc((size_t)(N*sizeof(float)));
	return m;
}

float* convert_matrix( float *a, long nrl, long nrh, long ncl, long nch)
{
	/* admpffpsjfpsdjpfsfpffssfpdofdpofdofsop */

    long i,j,k,N,
        nrow=nrh-nrl+1, ncol=nch-ncl+1;
    N = (nrow+1)*(ncol+1);
    float *m;

    m=(float *) malloc((size_t)(N*sizeof(float)));

    for (i = 1,k=0; i <= nrow; ++i)
        for (j = 1; j <= ncol; ++j,k++)
        	m[i*(ncol)+j]=a[k];

    if (!m) nrerror("allocation failure in convert_matrix()");
	 
    // for (i = 1; i <= nrow; ++i){
    //     printf("\n");
    //     for (j = 1; j <= ncol; ++j)
    //         printf("[%f]",m[i*ncol+j]);
    // }      
    return m;
}
