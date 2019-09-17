#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include"auxfunc.h"

void max_elem(int *piv_elem,int n,float *mat){
	int r,c;
	int max_i = 1;											//first coordenate i
	int max_j = 2;											//first coordenate j

	#pragma acc loop 
	for (r = 1; r <= n-1; r++) 
      	for (c = r+1; c <= n; c++)
      		if(fabs(mat[r*n+c]) > fabs(mat[max_i*n+max_j])){ //if exists a higher element
      			max_i = r;									//replace new coor
      			max_j = c;
    		}
    piv_elem[0] = max_i;									//store new coordenates
    piv_elem[1] = max_j;

}

float cal_tan(int max_i,int max_j,float *mat, int n){ 
	float num;
	float den;
	float a1;
	float a2;
	float a3;

	num = 2 * (mat[max_i*n+max_j]);								
	if(mat[max_i*n+max_i] < mat[max_i*n+max_i])
		num = -num;

	a1 = mat[max_i*n+max_i] - mat[max_j*n+max_j]; 
	a2 = a1*a1;
	a3 = 4 * mat[max_i*n+max_j]*mat[max_i*n+max_j];
	den = a2 + a3;
	den = sqrt(den);
	den = abs(a1) + den;
	return num/den;
}

float cal_cos(float tang){						//cos = 1/√(1+tan^2)
	float cose;
	cose = 1 + (tang * tang);
	cose = sqrt(cose);
	cose = 1 / cose;
	return cose;
}

float cal_sin(float cose, float tang){			//sin = cos*tan
	float sino;
	sino = cose*tang;
	return sino;
}

void mat_mult(int n,float *mat, float *T,float *mat_temp){
	int i,j,k;

	#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				//Premultiplication
    	for (j = 1 ; j <= n ; j++ ){
      		mat_temp[i*n+j] = 0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += T[k*n+i] * mat[k*n+j];
			}
		}
	}
	#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){					//Postmultiplication
    	for (j = 1 ; j <= n ; j++ ){
      		mat[i*n+j] = 0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat[i*n+j] += mat_temp[i*n+k] * T[k*n+j];
			}
		}
	}
}

void mat_mult2(int n,float *A, float *B,float *C){
	int i,j,k;
	for (i = 1 ; i <= n ; i++ ){ 			
    	for (j = 1 ; j <= n ; j++ ){
      		C[i*n+j] = 0.0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            C[i*n+j] += A[i*n+k] * B[k*n+j];
			}
		}
	}

}

void mat_mult_inv(int n,float *A, float *B,float *C){
	int i,j,k;
	for (i = 1 ; i <= n ; i++ ){ 			
    	for (j = 1 ; j <= n ; j++ ){
      		C[i*n+j] = 0.0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            C[i*n+j] += A[k*n+i] * B[k*n+j];
			}
		}
	}

}

void mult_eigenvec(int n,float *T,float *eigvec,float *mat_temp){
	int i,j,k;

	#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 								//eigenvec = eigenvec * T
    	for (j = 1 ; j <= n ; j++ ){
    		mat_temp[i*n+j] = 0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += eigvec[k*n+i] * T[k*n+j];
			}
		}
	}

	#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				
    	for (j = 1 ; j <= n ; j++ ){
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            eigvec[i*n+j] = mat_temp[i*n+j] ;
			}
		}
	}
}

void new_T_mat(int max_i, int max_j,int n,float *mat,float *T,float *mat_temp,float *eigvec){
	float tang, cose, sino;
	int c,r;

	tang = cal_tan(max_i,max_j,mat,n);
	cose = cal_cos(tang);
	sino = cal_sin(cose,tang);

	for (r = 1; r <= n; r++){				//Generate identity matrix
      	for (c = 1; c <= n; c++) 
     		T[r*n+c] = 0.0;		
      	T[r*n+r] = 1.0;	
	}
											//T Rotating matrix
    T[max_i*n+max_i] = cose;				
    T[max_j*n+max_j] = cose;
    T[max_i*n+max_j] = -sino; 				//Element to eliminate	
    T[max_j*n+max_i] = sino;		
    
}

void jacobiMultip (float *mat, int n, int ndm, float *eigvec, float eigval[],  int *nrot){

/*****************************+*******************************/	
	//On input
	//a: Contains the matrix to be diagonalized.
	//n: Order of matrix a.
	//ndim: Dimension of a.
	//eigvec: eigenvectors to be computed v

	//On output
	//eigval: Contains the eigenvalues in ascending order d
	//u: Contains the corresponding eigenvectors.
	//nrot: Number of Jacobi rotations.
/*****************************+*******************************/	

	int i,j;
	int *piv_elem;					//Keep coordenates of an elemnt i,j
	bool min = false;
	float EPS = .0000001;
	float *T;
	float *mat_temp; 

	mat_temp = (float *) vector(1,n*n);					
	T=(float *) vector(1,n*n);
	piv_elem=(int *) malloc (2 * sizeof (int));

	for (i = 1; i <= n; i++){						//Initializing Identity matrix
      	for (j = 1; j <= n; j++) 
     		eigvec[i*n+j] = 0.0;		
      	eigvec[i*n+i] = 1.0;	
	}
	
	for (*nrot = 0; min == false ; ++*nrot){
		max_elem(piv_elem,n,mat);	//Search for max element in tringle up

		if(fabs(mat[piv_elem[0]*n+piv_elem[1]]) < EPS || *nrot >= 100 ) //if max element doesnt exist more
			min=true;	
		
		else{
			new_T_mat(piv_elem[0],piv_elem[1],n,mat,T,mat_temp,eigvec); //Calculate T matrix
			mult_eigenvec(n,T,eigvec,mat_temp);							//Compute eigenvec
			//mat_mult(n,mat,T,mat_temp);									//Pre and Postmultiplicate T
			mat_mult_inv(n,T,mat,mat_temp);
			mat_mult2(n,mat_temp,T,mat);							
		}

		for (i = 1; i <= n; ++i)
			eigval[i]=mat[i*n+i];

	}
 
 	//*nrot = nrota;
 	//printf("rooooooot %d\n", *nrot);
	// free_vector(mat_temp,1,n*n);
	// free_vector(T,1,n);
	// free_vector(piv_elem,1,1);

}