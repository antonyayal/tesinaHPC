/*
 * jacobiACUDA.h
 * 
 *
 * José Antonio Ayala Barbosa
 * Oct, 2018
 *
 * 
 */

#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "auxFuncs.h"

__global__ void kernel_max_elem(int *piv_elem,int n,float *mat){
	int r,c;
	int max_i = 0;											//first coordenate i
	int max_j = 1;											//first coordenate j
	
	for (r = 0; r < n-1; r++) 
      	for (c = r+1; c < n; c++)
      		if(fabs(mat[r*n+c]) > fabs(mat[max_i*n+max_j])){ //if exists a higher element
      			max_i = r;									//replace new coor
      			max_j = c;
    		}
    piv_elem[0] = max_i;									//store new coordenates
    piv_elem[1] = max_j;

}

__global__ void kernel_set_max_elem(float *mat,int n, int *piv_elem,float *max_elem){
	*max_elem=mat[piv_elem[0]*n+piv_elem[1]];
}

__device__ float cal_tan(int max_i,int max_j,float *mat, int n){ 
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

__device__ float cal_cos(float tang){						//cos = 1/√(1+tan^2)
	float cose;
	cose = 1 + (tang * tang);
	cose = sqrt(cose);
	cose = 1 / cose;
	return cose;
}

__device__ float cal_sin(float cose, float tang){			//sin = cos*tan
	float sino;
	sino = cose*tang;
	return sino;
}

__global__ void kernel_new_T_mat(int *piv_elem,int n,float *mat,float *T,float *mat_temp){
	float tang, cose, sino;
	int c,r;
	int max_i = piv_elem[0];
	int max_j = piv_elem[1];

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if(row >= n|| col >= n) 
		return;

	for (row = row; row < n; ++row)
		T[row*n+col] = 0.0;
	T[col*n+col]=1.0;

	if(max_j==col){
		tang = cal_tan(max_i,max_j,mat,n);
		cose = cal_cos(tang);
		sino = cal_sin(cose,tang);

	    T[max_i*n+max_i] = cose;				
	    T[max_j*n+max_j] = cose;
	    T[max_i*n+max_j] = -sino; 				//Element to eliminate	
	    T[max_j*n+max_i] = sino;
	}	
}

__global__ void kernel_mat_mult(int n,float *A, float *B,float *C){

	float Cvalue;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	if(row >= n|| col >= n) 
		return;

	for (row = row; row < n; ++row){
		Cvalue = 0.0;
		for (e = 0; e < n; e++)
			Cvalue += (A[row * n + e]) * (B[e * n + col]);
		C[row * n + col] = Cvalue;
	}
}

__global__ void kernel_mat_mult_tra(int n,float *A, float *B,float *C){

	float Cvalue;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e,i;

	if(row >= n|| col >= n) 
		return;

	for (i = row; i < n; i++){
		Cvalue = 0.0;
		for (e = 0; e < n; e++)
			Cvalue += (A[e * n + i]) * (B[e * n + col]);
		C[i * n + col] = Cvalue;
	}
}

__global__ void kernel_copy_mat(int n,float *A, float *B){
	
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	if(row >= n|| col >= n) 
		return;

    for (e = 0; e < n; e++)
        B[col*n+e]=A[col*n+e];
}

__global__ void kernel_copy_diag(int n,float *A, float *B){
	
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if(row >= n|| col >= n) 
		return;

	B[col]=A[col*n+col];
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
	bool min = false;
	float EPS = .0000001;
	float *max_elem;
	float max_e;

	float *d_mat; 
	float *d_eigvec; 
	float *d_eigval; 
	float *d_T; 
	float *d_mat_temp;
	int *d_piv_elem;
	float *d_max_elem;

	size_t size = (n+1) * (n+1) * sizeof(float);
	dim3 dimGrid( 32 );         // 32 x 1 x 1
	dim3 dimBlock( 64 ); 		// 64 x 1 x 1 

	cudaMalloc(&d_mat, size);
	cudaMalloc(&d_eigvec, size);
	cudaMalloc(&d_eigval, n*sizeof (int));
	cudaMalloc(&d_T, size);
	cudaMalloc(&d_mat_temp, size);
	cudaMalloc(&d_piv_elem, 2*sizeof (int));
	cudaMalloc(&d_max_elem,sizeof (float));

	cudaMemcpy(d_mat, mat, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_eigvec, eigvec, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_eigval, eigval, n*sizeof (int), cudaMemcpyHostToDevice);
	
	for (*nrot = 0; min == false ; ++*nrot){

		kernel_max_elem<<<1,1>>>(d_piv_elem,n,d_mat);	//Search for max element in upper tringle 
		kernel_set_max_elem<<<1,1>>>(d_mat,n,d_piv_elem,d_max_elem);
		cudaMemcpy(&max_e, d_max_elem, sizeof(float), cudaMemcpyDeviceToHost);

		if(fabs(max_e) < EPS || *nrot >= 10000 ) //if max element doesnt exist more
			min=true;	
		
		else{
			kernel_new_T_mat<<<dimGrid, dimBlock>>>(d_piv_elem,n,d_mat,d_T,d_mat_temp);
			kernel_mat_mult<<<dimGrid, dimBlock>>>(n,d_eigvec,d_T,d_mat_temp);
			kernel_copy_mat<<<dimGrid, dimBlock>>>(n,d_mat_temp,d_eigvec);
			kernel_mat_mult_tra<<<dimGrid, dimBlock>>>(n,d_T,d_mat,d_mat_temp);
			kernel_mat_mult<<<dimGrid, dimBlock>>>(n,d_mat_temp,d_T,d_mat);
		}
	}

	kernel_copy_diag<<<dimGrid, dimBlock>>>(n,d_mat,d_eigval);

	cudaMemcpy(mat, d_mat, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(eigvec, d_eigvec, size, cudaMemcpyDeviceToHost); 
	cudaMemcpy(eigval,d_eigval,n*sizeof(float),cudaMemcpyDeviceToHost);

	cudaFree(d_mat); 
	cudaFree(d_eigvec); 
	cudaFree(d_eigval); 
	cudaFree(d_T); 
	cudaFree(d_mat_temp);
	cudaFree(d_piv_elem);
	cudaFree(d_max_elem);
}