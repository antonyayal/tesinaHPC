/* Driver for routine EIGSRT  */
#include <stdio.h>
//#include "nr.h"
//#include "nrutil.h"
#include <sys/time.h>
#include <stdlib.h>
//#include"jacobiMultipCUDA.h"

//#define NP 128



#include <stdbool.h>
#define NR_END 1
#define FREE_ARG char*

void get_walltime_(double* wcTime) {
   struct timeval tp;
   gettimeofday(&tp, NULL);
   *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void get_walltime(double* wcTime) {
   get_walltime_(wcTime);
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

float *matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long N,
        nrow=nrh-nrl+1, 
        ncol=nch-ncl+1;
    N = (nrow+1)*(ncol+1);
    float *m;

    m=(float *) malloc((size_t)(N*sizeof(float)));
	if (!m) nrerror("allocation failure in matrix()");
    return m;
}

float* convert_matrix( float *a, long nrl, long nrh, long ncl, long nch)
{
	/* admpffpsjfpsdjpfsfpffssfpdofdpofdofsop */

    long i,j,k,N,
        nrow=nrh-nrl+1, ncol=nch-ncl+1;
    N = (nrow+1)*(ncol+1);
    float *m;

    m=(float *) malloc((size_t)(N*sizeof(float)));

    for (i = 1,k=0; i <= nrow; ++i)
        for (j = 1; j <= ncol; ++j,k++)
        	m[i*(ncol)+j]=a[k];

    if (!m) nrerror("allocation failure in convert_matrix()");
	 
    // for (i = 1; i <= nrow; ++i){
    //     printf("\n");
    //     for (j = 1; j <= ncol; ++j)
    //         printf("[%f]",m[i*ncol+j]);
    // }      
    return m;
    
}

void free_vector(float *v,long nl,long nh){
/* free a float vector allocated with vector() */

    free((FREE_ARG) (v+nl-NR_END));
}

//void jacobiMultip (float *mat, int n, int ndm, float *eigvec, float eigval[],  int *nrot);







void max_elem(int *piv_elem,int n,float *mat){
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

	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				//Premultiplication
    	for (j = 1 ; j <= n ; j++ ){
      		mat_temp[i*n+j] = 0;
        	//#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += T[k*n+i] * mat[k*n+j];
			}
		}
	}

	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){					//Postmultiplication
    	for (j = 1 ; j <= n ; j++ ){
      		mat[i*n+j] = 0;
        	//#pragma acc loop
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
        	//#pragma acc loop
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
        	//#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            C[i*n+j] += A[k*n+i] * B[k*n+j];
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
        	//#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += eigvec[k*n+i] * T[k*n+j];
			}
		}
	}

	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				
    	for (j = 1 ; j <= n ; j++ ){
        	//#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            eigvec[i*n+j] = mat_temp[i*n+j] ;
			}
		}
	}
}

void new_T_mat(int max_i, int max_j,int n,float *mat,float *T){
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

__global__ void kernel_mat_mult_inv(int *n,float *A, float *B,float *C){

	float Cvalue = 0.0;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	//printf("n: %d\n", *n);

	if(row >= (*n) || col >= (*n)) 
		return;

	for (e = 1; e <= (*n); ++e)
		Cvalue += (A[e * (*n) + col]) * (B[e * (*n) + col]);
	
	C[row * (*n) + col] = Cvalue;
}

__global__ void kernel_mat_mult(int *n,float *A, float *B,float *C){

	float Cvalue = 0.0;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int e;

	//printf("n: %d\n", *n);

	if(row >= (*n) || col >= (*n)) 
		return;

	for (e = 1; e <= (*n); e++)
		Cvalue += (A[row * (*n) + e]) * (B[e * (*n) + col]);
	
	C[row * (*n) + col] = Cvalue;
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

	int *d_n;
	float *d_mat; 
	float *d_T; 
	float *d_mat_temp; 

	size_t size = (n+1) * (n+1) * sizeof(float);

	cudaError_t err;

	dim3 dimGrid( 512 );         // 512 x 1 x 1
	dim3 dimBlock( 1024, 1024 ); // 1024 x 1024 x 1 

	cudaMalloc(&d_n,sizeof (int));
	cudaMalloc(&d_mat, size);
	cudaMalloc(&d_T, size);
	cudaMalloc(&d_mat_temp, size);

	cudaMemcpy(d_n, &n, sizeof (int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mat_temp, &mat_temp, size, cudaMemcpyHostToDevice);

	//mat_temp = (float *) vector(1,n*n);					
	//T=(float *) vector(1,n*n);
	piv_elem=(int *) malloc (2 * sizeof (int));

	for (i = 1; i <= n; i++){						//Initializing Identity matrix
      	for (j = 1; j <= n; j++) 
     		eigvec[i*n+j] = 0.0;		
      	eigvec[i*n+i] = 1.0;	
	}
	//CPU
	for (*nrot = 0; min == false ; ++*nrot){
		max_elem(piv_elem,n,mat);	//Search for max element in tringle up

		if(fabs(mat[piv_elem[0]*n+piv_elem[1]]) < EPS || *nrot >= 100 ) //if max element doesnt exist more
			min=true;	
		
		else{
			new_T_mat(piv_elem[0],piv_elem[1],n,mat,T); //Calculate T matrix
			//mult_eigenvec(n,T,eigvec,mat_temp);							//Compute eigenvec
			mat_mult_inv(n,T,mat,mat_temp);
			mat_mult2(n,mat_temp,T,mat);
			
			// cudaMemcpy(d_mat, &mat, size, cudaMemcpyHostToDevice);
			// cudaMemcpy(d_T, &T, size, cudaMemcpyHostToDevice);
			// kernel_mat_mult_inv<<<dimGrid, dimBlock>>>(d_n,d_T,d_mat,d_mat_temp);
			// kernel_mat_mult<<<dimGrid, dimBlock>>>(d_n,d_mat_temp,d_T,d_mat);
			
			// cudaMemcpy(&mat, d_mat, size, cudaMemcpyDeviceToHost);

			
			//err=cudaMemcpy(&mat, d_mat, size, cudaMemcpyDeviceToHost);	
			//printf("Copy MAT off of device: %s\n",cudaGetErrorString(err));

			printf("\nRotación: %d\n",*nrot );
			for (i = 1; i <= n; ++i){
            printf("\n");
            for (j = 1; j <= n; ++j)
                printf("[%f]",mat[i*n+j]);
        }

		}

		for (i = 1; i <= n; ++i)
			eigval[i]=mat[i*n+i];

	}
 	cudaFree(d_mat);
  	cudaFree(d_T);
  	cudaFree(d_mat_temp);

  	//free(mat_temp);

 	//*nrot = nrota;
 	//printf("rooooooot %d\n", *nrot);
	// free_vector(mat_temp,1,n*n);
	// free_vector(T,1,n);
	// free_vector(piv_elem,1,1);

}




int main(int argc, char **argv)
{
    int NP;
    char *nombreArchivo=argv[1];
    double S,E;
	int i, j, nrot=0;
    FILE *archivo;
    float *c;

    if (fopen(nombreArchivo, "r") == NULL){
        printf("File not found\n");
        return 1;
    }else{
        archivo = fopen(nombreArchivo, "r");
        fscanf(archivo, "%d", &NP);
        c =(float *)matrix(1,NP-1,1,NP-1);
        for (i = 0; i < NP; i++){
            for (j = 0; j < NP; j++){
                fscanf(archivo, "%f", &c[i*NP+j]);
            }
        }
        fclose(archivo);
    }
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
              
        jacobiMultip(e,NP,NP,v,d,&nrot);
        
        get_walltime(&E);

        for (i = 1; i <= NP; ++i){
            printf("\n");
            for (j = 1; j <= NP; ++j)
                printf("[%f]",v[i*NP+j]);
        }
        
        printf("\nd\n");
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

    printf("Rotations: %d\n",nrot );

    printf("Total time:%f sec\n", E-S);
	

    //free_vector(d,1,NP);
    //free_vector(v,1,NP*NP);
    //free_vector(e,1,NP*NP);
    // free(d);
    // free(v);
    // free(e);

    return 0;
	
}
