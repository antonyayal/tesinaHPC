#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

__global__ void kernel_mat_mult(int *n,float *A, float *B,float *C){

	float Cvalue;
	int row = blockIdx.y * blockDim.y + threadIdx.y+1;
	int col = blockIdx.x * blockDim.x + threadIdx.x+1;
	int e;

	if(row > (*n)|| col > (*n)) 
		return;

	for (row = row; row < (*n); ++row){
		Cvalue = 0.0;
		for (e = 1; e < (*n); e++)
			Cvalue += (A[row * (*n) + e]) * (B[e * (*n) + col]);
		C[row * (*n) + col] = Cvalue;
	}
}

__global__ void kernel_mat_mult_tras(int *n,float *A, float *B,float *C){

	float Cvalue;
	int row = blockIdx.y * blockDim.y + threadIdx.y+1;
	int col = blockIdx.x * blockDim.x + threadIdx.x+1;
	int e,i;

	if(row > (*n)|| col > (*n)) 
		return;

	for (i = col; i < (*n); ++i){
		Cvalue = 0.0;
		for (e = 1; e < (*n); e++)
			Cvalue += (A[e * (*n) + i]) * (B[e * (*n) + col]);
		C[i * (*n) + col] = Cvalue;
	}

	// for (col = col; col < (*n); ++col,++row){
	// 	Cvalue = 0.0;
	// 	for (e = 1; e < (*n); e++)
	// 		Cvalue += (A[e * (*n) + col]) * (B[e * (*n) + col]);
	// 	C[row * (*n) + col] = Cvalue;
	// }
}

__global__ void imprimirMatDev(int *n,float *A,char ch){
	int i,j;
	printf("\nMatrix %c\n",ch);

	for (i = 1; i < *n; ++i){
            printf("\n");
            for (j = 1; j < *n; ++j)
                printf("[%f]",A[i*(*n)+j]);
    }
}

int main(int argc, char const *argv[])
{
	int n=4;
	int *d_n;
	float *d_mat_A;
	float *d_mat_B;
	float *d_mat_C;
	size_t size = (n+1) * (n+1) * sizeof(float);
	dim3 dimGrid( 32 );         // 512 x 1 x 1
	dim3 dimBlock( 64); // 1024 x 1024 x 1 

	float mat_A[4*4]= 
        {0.0,0.0,0.0,0.0,
         0.0,2.0,-1.0,0.0,
        0.0,-1.0,2.0,-1.0,
         0.0,0.0,-1.0,2.0,};
    float mat_B[4*4]= 
        {0.0,0.0,0.0,0.0,
         0.0,1.0,0.0,0.0,
         0.0,0.0,1.0,0.0,
         0.0,0.0,0.0,1.0,};
float mat_C[4*4]= 
        {0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,};
	cudaMalloc(&d_n,sizeof (int));
	cudaMalloc(&d_mat_A, size);
	cudaMalloc(&d_mat_B, size);
	cudaMalloc(&d_mat_C, size);

	printf("[%f]\n",mat_A[5] );
	printf("[%f]\n",mat_A[1*3+1] );

	cudaMemcpy(d_n, &n, sizeof (int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mat_A, mat_A, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mat_B, mat_B, size, cudaMemcpyHostToDevice);
	imprimirMatDev<<<1,1>>>(d_n,d_mat_A,'A');
	printf("\n");
	imprimirMatDev<<<1,1>>>(d_n,d_mat_B,'B');
	printf("\n");
	kernel_mat_mult<<<dimGrid, dimBlock>>>(d_n,d_mat_A,d_mat_B,d_mat_C);
	imprimirMatDev<<<1,1>>>(d_n,d_mat_C,'C');
	printf("\n");
	kernel_mat_mult_tras<<<dimGrid, dimBlock>>>(d_n,d_mat_A,d_mat_B,d_mat_C);
	imprimirMatDev<<<1,1>>>(d_n,d_mat_C,'D');
	cudaMemcpy(mat_C, d_mat_C, size, cudaMemcpyDeviceToHost);

	return 0;
}