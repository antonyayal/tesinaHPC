/* Driver for routine EIGSRT  */
#include <stdio.h>
//#include "nr.h"
//#include "nrutil.h"
#include <sys/time.h>
#include <stdlib.h>

//#define NP 128

int main(int argc, char **argv)
{
    int NP;
    char *nombreArchivo=argv[1];
    double S,E;
	int i, j, nrot=0;
    FILE *archivo;
    float *c;

    if (fopen(nombreArchivo, "r") == NULL){
        printf("File not found\n");
        return 1;
    }else{
        archivo = fopen(nombreArchivo, "r");
        fscanf(archivo, "%d", &NP);
        c =(float *)matrix(1,NP-1,1,NP-1);
        for (i = 0; i < NP; i++){
            for (j = 0; j < NP; j++){
                fscanf(archivo, "%f", &c[i*NP+j]);
            }
        }
        fclose(archivo);
    }
        float *d, *v, *e;
        
        d=(float *)vector(1,NP);
        v=(float *)matrix(1,NP,1,NP);
        e=(float *)convert_matrix(c,1,NP,1,NP);
        
        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%f]",e[i*NP+j]);
        }

        printf("\n****** Finding Eigenvectors ******\n");
        //jacobi(e,NP,d,v,&nrot);
        get_walltime(&S);
              
        jacobiMultip(e,NP,NP,v,d,&nrot);
        
        get_walltime(&E);

        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%f]",v[i*NP+j]);
        }
        
        printf("\nd\n");
        for (i = 1; i <= NP; ++i)
        {
            printf("[%f]",d[i]);
        }

        printf("\n******  Eigenvalues & Eigenvectors ******\n");
        for (i=1;i<=NP; i++) {
            printf("eigenvalue %3d, = %12.6f\n",i,d[i]);
            printf("eigenvector:\n");
            for (j=1;j<=NP; j++) {
                printf("%12.6f",v[i*NP+j]);
                if ((j % 5) == 0) printf("\n");
            }
                printf("\n");
        }

    printf("Rotations: %d\n",nrot );

    printf("Total time:%f sec\n", E-S);
	

    //free_vector(d,1,NP);
    //free_vector(v,1,NP*NP);
    //free_vector(e,1,NP*NP);

    return 0;
	
}
