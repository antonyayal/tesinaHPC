/*
 * driverA.cu
 * Driver for function JACOBI in CUDA.
 *
 * José Antonio Ayala Barbosa
 * Oct, 2018
 *
 * 
 */

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "auxFuncs.h"
#include "jacobiACUDA.h"

int main(int argc, char **argv)
{

	int n;
    char *nombreArchivo=argv[1];
    double S,E;
	int i, j, k, nrot=0;
    FILE *archivo;
    float *c,*mat,*eval,*evec;

    if (fopen(nombreArchivo, "r") == NULL){
        printf("File not found\n");
        return 1;
    }else{
        archivo = fopen(nombreArchivo, "r");
        fscanf(archivo, "%d", &n);
        c = vectorF(n*n);
        for (i = 0; i < n; ++i){
            for (j = 0; j < n; ++j){
                fscanf(archivo, "%f", &c[i*n+j]);
            }
        }
        fclose(archivo);
    }

    mat=vectorF(n*n);
    evec=vectorF(n*n);
    eval=vectorF(n);

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            mat[i*n+j]=c[i*n+j];

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++) 
            evec[i*n+j] = 0.0;        
        evec[i*n+i] = 1.0;    
    }

    for (i = 0; i < n; ++i){
            printf("\n");
            for (j = 0; j < n; ++j)
                printf("[%f]",c[i*n+j]);
    }

    printf("\n****** Finding Eigenvalues ******\n");
    	get_walltime(&S);
              
        jacobiMultip(mat,n,evec,eval,&nrot);
        
        get_walltime(&E);

    printf("\n******  Eigenvalues ******\n");
    for (i = 0; i < n; ++i){
            printf("\n");
            for (j = 0; j < n; ++j)
                printf("[%f]",mat[i*n+j]);
    }
    printf("\n");
	    for (i = 0; i < n; ++i){
	    	printf("eigenvalue %3d, = %12.6f\n",i+1,eval[i]);
            printf("eigenvector:\n");
            for (j=0,k=1;j<n; j++,k++) {
                printf("%12.6f",evec[i*n+j]);
                if ((k%5) == 0) printf("\n");
            }
            printf("\n");
        }

    printf("Rotations: %d\n",nrot );

    printf("Total time:%f sec\n", E-S);

    free(c);
    free(mat);
    free(evec);
    free(eval);
    
return 0;
}