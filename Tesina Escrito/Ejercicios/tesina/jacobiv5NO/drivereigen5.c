/* Driver for routine EIGSRT  */
#include <stdio.h>
#include "nr.h"
#include "nrutil.h"
#include <sys/time.h>
#include <stdlib.h>

#define NP 3

main()
{
    double S,E;
	int i, j, nrot;
    static float c[NP][NP]=
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
        
        float *d, **v, **e;
        float *le;
        float *lv;
        
        d=vector(1,NP);
        v=matrix(1,NP,1,NP);
        e=convert_matrix(&c[0][0],1,NP,1,NP);

        le=e;
        lv=(float *)v;

        printf("\n****** Original Matrix ******\n");
        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%12.6f]",le[i*NP+j]);
        }

        printf("\n****** Finding Eigenvectors ******\n");
        //jacobi(e,NP,d,v,&nrot);
        get_walltime(&S);
        jacobil(le,NP,NP,lv,d,&nrot);
        get_walltime(&E);

        printf("\n****** Procesed Matrix ******\n");
        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%12.6f]",v[i][j]);
        }
        printf("\n");
        for (i = 1; i <= NP; ++i)
        {
            printf("[%12.6f]",d[i]);
        }

        printf("\n******  Eigenvalues & Eigenvectors ******\n");
        for (i=1;i<=NP; i++) {
            printf("eigenvalue %3d, = %12.6f\n",i,d[i]);
            printf("eigenvector:\n");
            for (j=1;j<=NP; j++) {
                printf("%12.6f",v[i][j]);
                if ((j % 5) == 0) printf("\n");
            }
                printf("\n");
        }

    printf("Total time:%f sec\n", E-S);
	

    free_convert_matrix(e,1,NP,1,NP);
    free_matrix(v,1,NP,1,NP);
    free_vector(d,1,NP);
	
}
