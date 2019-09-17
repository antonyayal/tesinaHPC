/* Driver for routine EIGSRT  */
#include <stdio.h>
//#include "nr.h"
//#include "nrutil.h"
#include <sys/time.h>
#include <stdlib.h>

#define NP 3

main()
{
    double S,E;
	int i, j, nrot;
    static float c[NP*NP]=
        {2.0,-1.0,0.0,
         -1.0,2.0,-1.0,
         0.0,-1.0,2.0,};    
    // static float c[5][5]=
    // {-2.0,-1.0,0.0,1.0,2.0,
    //  -1.0,-1.0,0.0,1.0,2.0,
    //  0.0,0.0,0.0,1.0,2.0,
    //  1.0,1.0,1.0,1.0,2.0,
    //  2.0,2.0,2.0,2.0,2.0,};
    // static float c [NP][NP]=
    // {5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,-4.0,
    //  4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,
    //  3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,
    //  2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,
    //  1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,
    //  0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,
    //  -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,
    //  -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,
    //  -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,
    //  -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,};
        
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
        jacobil(e,NP,NP,v,d,&nrot);
        get_walltime(&E);

        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%f]",v[i*NP+j]);
        }
        printf("\n");
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

    printf("Total time:%f sec\n", E-S);
	

    free_vector(d,1,NP);
    free_vector(v,1,NP*NP);
    free_vector(e,1,NP*NP);
	
}
