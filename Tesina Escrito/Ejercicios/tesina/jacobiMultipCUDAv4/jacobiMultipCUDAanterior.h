#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//#include"auxfunc.h"

#define BLOCK_SIZE 16

__device__ float tang;
__device__ float cose;
__device__ float sino;

__global__ void print(){
	printf("Imprimiendo.........\n");
}
__global__ void kernel_max_elem(int *piv_elem,int *n,float *mat,float *m_e){
	int r,c;
	int max_i = 1;											//first coordenate i
	int max_j = 2;											//first coordenate j

	//#pragma acc loop 
	for (r = 1; r <= *n-1; r++) 
      	for (c = r+1; c <= *n; c++)
      		if(fabs(mat[r*(*n)+c]) > fabs(mat[max_i*(*n)+max_j])){ //if exists a higher element
      			max_i = r;									//replace new coor
      			max_j = c;
    		}
    piv_elem[0] = max_i;									//store new coordenates
    piv_elem[1] = max_j;

    *m_e = mat[max_i*(*n)+max_j];

    printf("Pasoooooo dentro %f\n",*m_e);
}

float max_elem(int *piv_elem,int n,float *mat){//,float *m_e){
	int r,c;
	int max_i = 1;											//first coordenate i
	int max_j = 2;											//first coordenate j

	//#pragma acc loop 
	for (r = 1; r <= n-1; r++) 
      	for (c = r+1; c <= n; c++)
      		if(fabs(mat[r*n+c]) > fabs(mat[max_i*n+max_j])){ //if exists a higher element
      			max_i = r;									//replace new coor
      			max_j = c;
    		}
    piv_elem[0] = max_i;									//store new coordenates
    piv_elem[1] = max_j;

    return mat[max_i*n+max_j];

    //printf("Pasoooooo dentro %f\n",*m_e);
}

__global__ void kernel_cal_tan(int max_i,int max_j,float *mat, int *n){ 
	float num;
	float den;
	float a1;
	float a2;
	float a3;
	//int num = *n;

	num = 2 * (mat[max_i*(*n)+max_j]);								
	if(mat[max_i*(*n)+max_i] < mat[max_i*(*n)+max_i])
		num = -num;

	a1 = mat[max_i*(*n)+max_i] - mat[max_j*(*n)+max_j]; 
	a2 = a1*a1;
	a3 = 4 * mat[max_i*(*n)+max_j]*mat[max_i*(*n)+max_j];
	den = a2 + a3;
	den = sqrt(den);
	den = abs(a1) + den;
	tang = num/den;
}

__global__ void kernel_cal_cos(){						//cos = 1/√(1+tan^2)
	float cose;
	cose = 1 + (tang * tang);
	cose = sqrt(cose);
	cose = 1 / cose;
	//cose;
}

__global__ void kernel_cal_sin(){			//sin = cos*tan
	//float sino;
	sino = cose*tang;
	//return sino;
}

__global__ void kernel_mat_mult(int *n,float *A, float *B,float *C){

	float Cvalue = 0.0;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	if(row >= (*n) || col >= (*n)) 
		return;

	for (e = 1; e <= (*n); ++e)
		Cvalue += (A[row * (*n) + e]) * (B[e * (*n) + col]);
	
	C[row * (*n) + col] = Cvalue;
}

__global__ void mat_mult(int *n,float *mat, float *T,float *mat_temp){
	int i,j,k;

	//#pragma acc loop
	for (i = 1 ; i <= (*n) ; i++ ){ 				//Premultiplication
    	for (j = 1 ; j <= (*n) ; j++ ){
      		mat_temp[i*(*n)+j] = 0;
        	//#pragma acc loop
        	for (k = 1 ; k <= (*n) ; k++ ){
	            mat_temp[i*(*n)+j] += T[k*(*n)+i] * mat[k*(*n)+j];
			}
		}
	}
	//#pragma acc loop
	for (i = 1 ; i <= (*n) ; i++ ){					//Postmultiplication
    	for (j = 1 ; j <= (*n) ; j++ ){
      		mat[i*(*n)+j] = 0;
        	//#pragma acc loop
        	for (k = 1 ; k <= (*n) ; k++ ){
	            mat[i*(*n)+j] += mat_temp[i*(*n)+k] * T[k*(*n)+j];
			}
		}
	}
}

void mult_eigenvec(int n,float *T,float *eigvec,float *mat_temp){
	int i,j,k;

	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 								//eigenvec = eigenvec * T
    	for (j = 1 ; j <= n ; j++ ){
    		mat_temp[i*n+j] = 0;
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += eigvec[k*n+i] * T[k*n+j];
			}

			//eigvec[i*n+j] = mat_temp[i*n+j] ;

		}
	}

	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				
    	for (j = 1 ; j <= n ; j++ ){
        	#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            eigvec[i*n+j] = mat_temp[i*n+j] ;
			}
		}
	}
}

__global__ void kernel_new_T_mat(int max_i, int max_j,int *n,float *mat,float *T,float *mat_temp,float *eigvec){
	//float tang, cose, sino;
	int c,r;

	for (r = 1; r <= (*n); r++){				//Generate identity matrix
      	for (c = 1; c <= (*n); c++) 
     		T[r*(*n)+c] = 0.0;		
      	T[r*(*n)+r] = 1.0;	
	}
											//T Rotating matrix
    T[max_i*(*n)+max_i] = cose;				
    T[max_j*(*n)+max_j] = cose;
    T[max_i*(*n)+max_j] = -sino; 				//Element to eliminate	
    T[max_j*(*n)+max_i] = sino;		
    
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
	// int nrota;
	int *piv_elem;					//Keep coordenates of an elemnt i,j
	//float sino,cose,tang;
	bool min = false;
	//bool max = true;
	float EPS = .0000001;
	float *T;
	float *mat_temp;
	float m_e;

	int *d_n;
	float *d_m_e;
	int *d_piv_elem;
	float *d_mat; 
	float *d_T; 
	float *d_mat_temp; 
	float *d_eigvec; 

	size_t sizeInt = sizeof (int);
	size_t sizeFloat = sizeof (float);

	size_t size = (n+1) * (n+1) * sizeof(float);

	dim3 dimGrid( 512 );         // 512 x 1 x 1
	dim3 dimBlock( 1024, 1024 ); // 1024 x 1024 x 1 

	mat_temp = (float *) vector(1,n*n);					
	T=(float *) vector(1,n*n);
	piv_elem=(int *) malloc (2 * sizeof (int));

	cudaMalloc(&d_piv_elem, 2 * sizeInt);
	cudaMalloc(&d_mat, size);
	cudaMalloc(&d_T, size);
	cudaMalloc(&d_mat_temp, size);
	cudaMalloc(&d_eigvec, size);

	for (i = 1; i <= n; i++){						//Initializing Identity matrix
      	for (j = 1; j <= n; j++) 
     		eigvec[i*n+j] = 0.0;		
      	eigvec[i*n+i] = 1.0;	
	}

	cudaMemcpy(d_n, &n, sizeInt, cudaMemcpyHostToDevice);
	//cudaMemcpy(d_m_e, &m_e, sizeFloat, cudaMemcpyHostToDevice);
	

	//cudaMemcpy(d_piv_elem, &piv_elem, 2*sizeInt, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mat, &mat, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_T, &T, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mat_temp, &mat_temp, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_eigvec, &eigvec, size, cudaMemcpyHostToDevice);
	
	for (*nrot = 0; min == false ; ++*nrot){
		//print<<<1,1>>>();
		// kernel_max_elem<<<1,1>>>(d_piv_elem,d_n,d_mat,d_m_e);	//Search for max element in tringle up
		// cudaMemcpy(&m_e, d_m_e, sizeFloat, cudaMemcpyDeviceToHost);
		cudaMemcpy(&mat, d_mat, size, cudaMemcpyDeviceToHost);
		m_e=max_elem(piv_elem,n,mat);
		cudaMemcpy(d_piv_elem, &piv_elem, sizeFloat, cudaMemcpyHostToDevice);
		printf("Pasoooooo %f\n",m_e);

		if(fabs(m_e) < EPS || *nrot >= 100 ) //if max element doesnt exist more
			min=true;	
		
		else{
			//new_T_mat(piv_elem[0],piv_elem[1],n,mat,T,mat_temp,eigvec); //Calculate T matrix
			//mult_eigenvec(n,T,eigvec,mat_temp);							//Compute eigenvec
			//mat_mult(n,mat,T,mat_temp);									//Pre and Postmultiplicate T
			
			kernel_cal_tan<<<1,1>>>(d_piv_elem[0],d_piv_elem[1],d_mat,d_n);
			kernel_cal_cos<<<1,1>>>();
			kernel_cal_sin<<<1,1>>>();
			kernel_new_T_mat<<<1,1>>>(d_piv_elem[0],d_piv_elem[1],d_n,d_mat,d_T,d_mat_temp,d_eigvec);
			mat_mult<<<1, 1>>>(d_n,d_mat,d_T,d_mat_temp);
			//kernel_mat_mult<<<dimGrid, dimBlock>>>(n,d_mat,d_T,d_mat_temp);
			//kernel_mat_mult<<<dimGrid, dimBlock>>>(n,d_mat_temp,d_T,d_mat);
										
		}

		//for (i = 1; i <= n; ++i)
			//eigval[i]=mat[i*n+i];

	}

	//cudaMemcpy(C.elements, d_C.elements, size, cudaMemcpyDeviceToHost);

	cudaMemcpy(&mat, d_mat, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(&T, d_T, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(&mat_temp, d_mat_temp, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(&eigvec, d_eigvec, size, cudaMemcpyDeviceToHost);
 	
 	cudaFree(d_mat);
  	cudaFree(d_T);
  	cudaFree(d_mat_temp);
  	cudaFree(d_eigvec);
  	// free(mat_temp);
  	// free(T);
  	// free(piv_elem);

	// free_vector(mat_temp,1,n*n);
	// free_vector(T,1,n);
	// free_vector(piv_elem,1,1);

}