/*
 * jacobiA.h
 * 
 *
 * José Antonio Ayala Barbosa
 * Oct, 2018
 *
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include"auxFuncs.h"

void max_elem(int *piv_elem,int n,float *mat){
	int r,c;
	int max_i = 0;											//first coordenate i
	int max_j = 1;											//first coordenate j

	#pragma acc parallel num_gangs(32), vector_length(64)
	{
		#pragma acc loop 
		for (r = 0; r < n-1; r++) 
	      	for (c = r+1; c < n; c++)
	      		if(fabs(mat[r*n+c]) > fabs(mat[max_i*n+max_j])){ //if exists a higher element
	      			max_i = r;									//replace new coor
	      			max_j = c;
	    		}
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

void new_T_mat(int max_i, int max_j,int n,float *mat,float *T,float *mat_temp){
	float tang, cose, sino;
	int c,r;

	tang = cal_tan(max_i,max_j,mat,n);
	cose = cal_cos(tang);
	sino = cal_sin(cose,tang);

	for (r = 0; r < n; r++){				//Generate identity matrix
      	for (c = 0; c < n; c++) 
     		T[r*n+c] = 0.0;		
      	T[r*n+r] = 1.0;	
	}
											//T Rotating matrix
    T[max_i*n+max_i] = cose;				
    T[max_j*n+max_j] = cose;
    T[max_i*n+max_j] = -sino; 				//Element to eliminate	
    T[max_j*n+max_i] = sino;		
    
}

void mat_mult(int n,float *A, float *B,float *C){
	int i,j,k;
    #pragma acc parallel num_gangs(32), vector_length(64)
	{
		#pragma acc loop
		for (i = 0 ; i < n ; i++ ){ 			
	    	for (j = 0 ; j < n ; j++ ){
	      		C[i*n+j] = 0.0;
	        	#pragma acc loop
	        	for (k = 0 ; k < n ; k++ ){
		            C[i*n+j] += A[i*n+k] * B[k*n+j];
				}
			}
		}
	}

}

void mat_mult_tra(int n,float *A, float *B,float *C){
	int i,j,k;
    #pragma acc parallel num_gangs(32), vector_length(64)
	{
   		#pragma acc loop
		for (i = 0 ; i < n ; i++ ){ 			
	    	for (j = 0 ; j < n ; j++ ){
	      		C[i*n+j] = 0.0;
	        	#pragma acc loop
	        	for (k = 0 ; k < n ; k++ ){
		            C[i*n+j] += A[k*n+i] * B[k*n+j];
				}
			}
		}
	}
}

void copy_mat(int n,float *A, float *B){
	int i,j;
	#pragma acc parallel num_gangs(32), vector_length(64)
	{
		#pragma acc loop
		for (i = 0; i < n; ++i)
	        for (j = 0; j < n; ++j)
	            B[i*n+j]=A[i*n+j];
    }
}


void jacobiMultip (float *mat, int n,float *eigvec, float *eigval,  int *nrot){

/*****************************+*******************************/	
	//On input
	//mat: Contains the matrix to be diagonalized.
	//n: Order of matrix a.

	//On output
	//eigvec: eigenvectors to be computed v
	//eigval: Contains the eigenvalues in ascending order
	//nrot: Number of Jacobi rotations.
/*****************************+*******************************/	

	int i,j;
	int *piv_elem;					//Keep coordenates of an elemnt i,j
	bool min = false;
	float EPS = .0000001;
	float *T;						//Contains the ratation matrix
	float *mat_temp; 				//A temporal matrix 

	mat_temp = vectorF(n*n);					
	T = vectorF(n*n);
	piv_elem = vectorI(2);
	
	for (*nrot = 0; min == false ; ++*nrot){
		max_elem(piv_elem,n,mat);	//Search for max element in tringle up

		if(fabs(mat[piv_elem[0]*n+piv_elem[1]]) < EPS || *nrot >= 50000 ) //if max element doesnt exist more
			min=true;	
		
		else{
			new_T_mat(piv_elem[0],piv_elem[1],n,mat,T,mat_temp); //Calculate T 
			mat_mult(n,eigvec,T,mat_temp);					//Eigenvect
			copy_mat(n,mat_temp,eigvec);
			mat_mult_tra(n,T,mat,mat_temp);
			mat_mult(n,mat_temp,T,mat);									
	}

		for (i = 0; i < n; ++i)
			eigval[i]=mat[i*n+i];
	}
}