#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <sys/time.h>

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void get_walltime_(double* wcTime) {
   struct timeval tp;
   gettimeofday(&tp, NULL);
   *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void get_walltime(double* wcTime) {
   get_walltime_(wcTime);
}

float *vectorF(int n){
	float *v;
	v=(float *) malloc((size_t)(n*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
    return v;
}

int *vectorI(int n){
	int *v;
	v=(int *) malloc((size_t)(n*sizeof(int)));
	if (!v) nrerror("allocation failure in vector()");
    return v;
}