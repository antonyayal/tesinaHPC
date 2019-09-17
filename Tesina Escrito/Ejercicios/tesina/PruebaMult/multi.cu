#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

__global__ void kernel_mat_mult(int *n,float *A, float *B,float *C){

	float Cvalue;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	if(row > (*n)|| col > (*n)) 
		return;

	for (row = row; row < (*n); ++row){
		Cvalue = 0.0;
		for (e = 0; e < (*n); e++)
			Cvalue += (A[row * (*n) + e]) * (B[e * (*n) + col]);
		C[row * (*n) + col] = Cvalue;
	}
}

__global__ void imprimirMatDev(int *n,float *A,char ch){
	int i,j;
	printf("\nMatrix %c\n",ch);
	for (i = 0; i < *n; ++i){
            printf("\n");
            for (j = 0; j < *n; ++j)
                printf("[%f]",A[i*(*n)+j]);
    }
}

int main(int argc, char const *argv[])
{
	int n=3;
	int *d_n;
	float *d_mat_A;
	float *d_mat_B;
	float *d_mat_C;
	size_t size = (n+1) * (n+1) * sizeof(float);
	dim3 dimGrid( 32 );         // 512 x 1 x 1
	dim3 dimBlock( 64); // 1024 x 1024 x 1 

	float mat_A[3*3]= 
        {2.0,-1.0,0.0,
        -1.0,2.0,-1.0,
         0.0,-1.0,2.0,};
    float mat_B[3*3]= 
        {2.0,0.0,0.0,
         0.0,2.0,0.0,
         0.0,0.0,2.0,};
float mat_C[3*3]= 
        {0.0,0.0,0.0,
         0.0,0.0,0.0,
         0.0,0.0,0.0,};
	cudaMalloc(&d_n,sizeof (int));
	cudaMalloc(&d_mat_A, size);
	cudaMalloc(&d_mat_B, size);
	cudaMalloc(&d_mat_C, size);


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
	cudaMemcpy(mat_C, d_mat_C, size, cudaMemcpyDeviceToHost);

	return 0;
}