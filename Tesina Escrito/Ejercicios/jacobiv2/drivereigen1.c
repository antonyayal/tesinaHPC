/* Driver for routine EIGSRT  */
#include <stdio.h>
//#include "nr.h"
//#include "nrutil.h"
#include <sys/time.h>

#define NP  5

void get_walltime_(double* wcTime) {
   struct timeval tp;
   gettimeofday(&tp, NULL);
   *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void get_walltime(double* wcTime) {
   get_walltime_(wcTime);
}
main()
{
   /* clock_t start = clock();*/
    double S,E;
	int i, j, nrot;
        static float c[3][3]=
            {2.0,-1.0,0.0,
             -1.0,2.0,-1.0,
             0,-1,2,};
        /*static float c[5][5]=
            {-2.0,-1.0,0.0,1.0,2.0,
            -1.0,-1.0,0.0,1.0,2.0,
            0.0,0.0,0.0,1.0,2.0,
            1.0,1.0,1.0,1.0,2.0,
            2.0,2.0,2.0,2.0,2.0,};*/
        /*static float c [NP][NP]=
	         {5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,-4.0,
             4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,
             3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,
             2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,
             1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,
             0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,
             -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,
             -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,
             -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,
             -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,};*/

        float *d, **v, **e;
        
        d=vector(1,NP);
        v=matrix(1,NP,1,NP);
        e=convert_matrix(&c[0][0],1,NP,1,NP);
        printf("****** Finding Eigenvectors ******\n");
        //jacobi(e,NP,d,v,&nrot);
        get_walltime(&S);
        jacobil(e,NP,NP,v,d,&nrot);
        get_walltime(&E);
        printf("****** unsorted eigenvectors ******\n");
        for (i=1;i<=NP; i++) {
    	    printf("eigenvalue %3d, = %12.6f\n",i,d[i]);
            printf("eigenvector:\n");
        	for (j=1;j<=NP; j++) {
        		printf("%12.6f",v[j][i]);
        		if ((j % 5) == 0) printf("\n");
        	}
                printf("\n");
	}

	printf("****** Sorting Eigenvectors ******\n\n");
	eigsrt(d,v,NP);
        
	printf("sorted eigenvectors:\n");
        for (i=1;i<=NP;i++){
		printf("eigenvalue %3d, = %12.6f\n",i,d[i]);
		printf("eigenvector:\n");
		for (j=1;j<=NP;j++)  {
			 printf("%12.6f",v[j][i]);
			  if ((j % 5) == 0) printf("\n");
		}
		printf("\n");
	}
    //printf("Total time: %f\n",((double)clock()-start)/CLOCKS_PER_SEC);
    printf("Total time:%f sec\n", E-S);
	free_convert_matrix(e,1,NP,1,NP);
	free_matrix(v,1,NP,1,NP);
	free_vector(d,1,NP);
}
