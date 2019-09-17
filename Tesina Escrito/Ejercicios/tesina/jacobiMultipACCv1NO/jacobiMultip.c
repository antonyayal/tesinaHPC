#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void max_elem(int *piv_elem,int n,float *mat,float EPS){
	int r,c;
	int max_i = 1;
	int max_j = 2;

	// for (r = 1; r <= n-1; r++) 
 //      	for (c = r+1; c <= n; c++)
 //      		if (fabs(mat[r*n+c])>EPS){
 //      			max_i = r;
	// 			max_j = c;
	// 			break;
 //      		}

 //    if (max_i==0)
 //    	return false;

	#pragma acc loop
	for (r = 1; r <= n-1; r++) 
		//#pragma acc loop
      	for (c = r+1; c <= n; c++)
      		if(fabs(mat[max_i*n+max_j])<fabs(mat[r*n+c])){
      			max_i = r;
      			max_j = c;
    		}
    piv_elem[0] = max_i;
    piv_elem[1] = max_j;

    //return true;
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

float cal_cos(float tang){
	float cose;
	cose = 1 + (tang * tang);
	cose = sqrt(cose);
	cose = 1 / cose;
	return cose;
}

float cal_sin(float cose, float tang){
	float sino;
	sino = cose*tang;
	return sino;
}

void mat_mult(int n,float *mat, float *T,float *mat_temp){
	int i,j,k;
	#pragma acc loop
	//#pragma acc data copy(mat), copy(T), copy(mat_temp)
	for (i = 1 ; i <= n ; i++ ){ 				//Premultiplication
    	//#pragma acc loop//kernels loop gang(100), vector(128)
    	//#pragma acc kernels
    	for (j = 1 ; j <= n ; j++ ){
      		mat_temp[i*n+j] = 0;
      		//#pragma acc kernels loop
      		#pragma acc loop //gangs(8) vector(32)
        	for (k = 1 ; k <= n ; k++ ){
	            mat_temp[i*n+j] += T[k*n+i] * mat[k*n+j];
			}
		}
	}
	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){					//Postmultiplication
    	//#pragma acc loop
    	for (j = 1 ; j <= n ; j++ ){
      		mat[i*n+j] = 0;
      		#pragma acc loop
        	for (k = 1 ; k <= n ; k++ ){
	            mat[i*n+j] += mat_temp[i*n+k] * T[k*n+j];
			}
		}
	}
}

void mult_eigenvec(int n,float *T,float *eigvec,float *mat_temp){
	int i,j,k;
	//#pragma acc loop
	for (i = 1 ; i <= n ; i++ ){ 				
		//#pragma acc loop
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
	//#pragma acc loop			
    	for (j = 1 ; j <= n ; j++ ){
    		//#pragma acc loop

        	for (k = 1 ; k <= n ; k++ ){
	            eigvec[i*n+j] = mat_temp[i*n+j] ;
			}
		}
	}
}

void new_T_mat(int max_i, int max_j,int n,float *mat,float *T,float *mat_temp,float *eigvec){
	float tang, cose, sino;
	int c,r,i,j;

	tang = cal_tan(max_i,max_j,mat,n);
	cose = cal_cos(tang);
	sino = cal_sin(cose,tang);

	for (r = 1; r <= n; r++){
      	for (c = 1; c <= n; c++) 
     		T[r*n+c] = 0.0;		
      	T[r*n+r] = 1.0;	
	}

    T[max_i*n+max_i] = cose;
    T[max_j*n+max_j] = cose;
    T[max_i*n+max_j] = -sino; //Termino a eliminar	
    T[max_j*n+max_i] = sino;		
    
  //   for (i = 1; i <= n; ++i){
  //       printf("\n");
  //       for (j = 1; j <= n; ++j)
  //           printf("[%f]",T[i*n+j]);
 	// }
    mult_eigenvec(n,T,eigvec,mat_temp);
	mat_mult(n,mat,T,mat_temp);
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
	int nrota;
	int *piv_elem;					//Keep coordenates of an elemnt i,j
	float sino,cose,tang;
	bool min = false;
	bool max = true;
	float EPS = .0000001;
	float *T;
	float *mat_temp; 

	mat_temp = (float *) vector(1,n*n);						//Temporal matrix for oper
	T=(float *) vector(1,n*n);								//Rotation matrix
	piv_elem=(int *) malloc (2 * sizeof (int));				//Coordenates of max element

	for (i = 1; i <= n; i++){								//Initializing eigenvector
      	for (j = 1; j <= n; j++) 
     		eigvec[i*n+j] = 0.0;		
      	eigvec[i*n+i] = 1.0;	
	}
	
	for (*nrot = 0; min == false ; *nrot++){				
		max_elem(piv_elem,n,mat,EPS);						//Find max element in tringle up

		if(fabs(mat[piv_elem[0]*n+piv_elem[1]]) < EPS || *nrot > 50 )	//If element almost 0
			min=true;	
		
		else{
		new_T_mat(piv_elem[0],piv_elem[1],n,mat,T,mat_temp,eigvec);		//
		max_elem(piv_elem,n,mat,EPS);
		}

		for (j = 1; j <= n; ++j)
			eigval[j]=mat[j*n+j];

		
		// if(fabs(mat[piv_elem[0]*n+piv_elem[1]]) < EPS || nrot > 50 )
		// 	min=true;	
	}
 
 	//*nrot = nrota;
 	//printf("rooooooot %d\n", *nrot);
	//free_vector(mat_temp,1,n*n);
	// free_vector(T,1,n);
	// free_vector(piv_elem,1,1);

}