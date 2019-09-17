#include<stdio.h>
#include<stdlib.h>
#include "auxfunc.h"


void mult(int n,float *A,float *B,float *C){
	int i,j,k;

	for (i = 1 ; i <= n ; i++ ){ 								//eigenvec = eigenvec * T
    	for (j = 1 ; j <= n ; j++ ){
    		C[i*n+j] = 0;
        	for (k = 1 ; k <= n ; k++ ){
	            C[i*n+j] += A[k*n+i] * B[k*n+j];
			}
		}
	}

}

int main(int argc, char const *argv[]){
	int i,j;
	int n = 3;

	float *A=(float *)vector(1,3);
	float *B=(float *)vector(1,3);
	float *C=(float *)vector(1,3);

	A={	0,0,0,0,
				0,1,0,0,
				0,0,1,0,
				0,0,0,1,};

	B={	0,0,0,0,
				0,1,2,3,
				0,4,5,6,
				0,7,8,9,};

	C={	0,0,0,0,
				0,1,0,0,
				0,0,1,0,
				0,0,0,1,};

	mult(n,A,B,C);

	for (i = 1 ; i <= n ; i++ ){ 
		printf("\n");
    	for (j = 1 ; j <= n ; j++ ){
    		printf("[%f]", C[i*n+j]);
    	}
    }

	return 0;
}