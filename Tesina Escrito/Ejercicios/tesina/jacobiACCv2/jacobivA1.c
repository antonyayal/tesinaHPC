/*************************************/
// jacobiv3.c Change matrix to major-row
//
/*************************************/
#include <math.h>

//#include <openacc.h>

#define NRANSI
//#include "nrutil.h"

#define ROTATE(a,n,i,j,k,l) \
		EPS=a[i*n+j];\
		temp=a[k*n+l];\
		a[i*n+j]=EPS-s*(temp+EPS*tau);\
		a[k*n+l]=temp+s*(EPS-temp*tau);

void jacobil(float *a, int n, int ndm, float *eigvec, float eigval[],  int *nrot){

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
	int j,iq,ip,i; 									//Counters.
	float *b,*z;									//Vectors
	float s,c,t;									//Trigonometric functions
	float theta,tau;								//
	float tresh,sm,temp,EPS;						//Auxiliar vars

	b=(float *) vector(1,n);
	z=(float *) vector(1,n);

		//#pragma acc loop
		for (ip=1;ip<=n;ip++) { 						//Initializating identity matrix.
			for (iq=1;iq<=n;iq++) 
				eigvec[ip*n+iq]=0.0;
			eigvec[ip*n+ip]=1.0;
	}
		//#pragma acc loop
		for (ip=1;ip<=n;ip++) {						
			b[ip]=eigval[ip]=a[ip*n+ip];				//Initializating b and eigvec to the diagonal af a.
			z[ip]=0.0;									//This vector will accumulate terms of the form tapq as in equation. (creo que old vector)
		}
	
	*nrot=0;	
	
							
	for (i=1;i<=50;i++) {
		//NO#pragma acc kernels	
		{
		sm=0.0;	
		//#pragma acc loop reduction(sm)				
		for (ip=1;ip<=n-1;ip++) {					//Sum off-diagonal upper triangle elements.
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip*n+iq]);
		}

		if (sm == 0.0) {
			eigsrt(eigval,eigvec,n);				//Sorting Eigenvalues & Eigenvectors
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}

		if (i < 4)									//First three steps.
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;	

		// printf("tresh %f  \n", tresh);
						
		#pragma acc loop
		for (ip=1;ip<=n-1;ip++) {

			for (iq=ip+1;iq<=n;iq++) {
				EPS=100.0*fabs(a[ip*n+iq]); 			//Epsilon of the machine
				// printf("EPS %f  \n", EPS);

				if (i > 4 							//After 4 sweeps, skip rotation if the off-diagonal element is too small. 
					&& (float)(fabs(eigval[ip])+EPS) == (float)fabs(eigval[ip])
					&& (float)(fabs(eigval[iq])+EPS) == (float)fabs(eigval[iq]))
					a[ip*n+iq]=0.0;					//As element is too small, = 0
				
				else if (fabs(a[ip*n+iq]) > tresh){ 	//If still can be computed.
					
					temp=eigval[iq]-eigval[ip];

					if ((float)(fabs(temp)+EPS) == (float)fabs(temp))		/******CHECAR******/
						t=(a[ip*n+iq])/temp;							//tanΘ = 1/(2Θ)
					
					else {
						theta=0.5*temp/(a[ip*n+iq]);					//Θ = (aqq-app)/2apq
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));	//tanΘ = ±1/(|Θ|+√(Θ^2+1))
						if (theta < 0.0) t = -t;					//plus sign is used if app≥aqq, else minus sign
					}

					c=1.0 / sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					temp=t*a[ip*n+iq];
					z[ip] -= temp;
					z[iq] += temp;
					eigval[ip] -= temp; 			//Row
					eigval[iq] += temp;				//Column
					a[ip*n+iq]=0.0;
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
					// printf("ROTATE 1\n");
					#pragma acc loop
					for (j=1;j<=ip-1;j++) {			//Sweep rows
						// printf("ip: %d\n",ip);
						// printf("iq: %d\n",iq);
						// printf("j: %d\n",j);
						// printf("a[%d][%d] a[%d][%d]\n",j,ip,j,iq);
						ROTATE(a,n,j,ip,j,iq)
					}
					// printf("ROTATE 2\n");
					#pragma acc loop
					for (j=ip+1;j<=iq-1;j++) {		//Sweep diagonals
						// printf("ip: %d\n",ip);
						// printf("iq: %d\n",iq);
						// printf("j: %d\n",j);
						// printf("a[%d][%d] a[%d][%d]\n",ip,j,j,iq);
						ROTATE(a,n,ip,j,j,iq)
						// for (int k = 1; k <= n; ++k){
				  //           printf("\n");
				  //           for (int l = 1; l <= n; ++l)
			   //              printf("[%f]  ",a[k*n+l]);
			        	
      //   				}
					}
					// printf("ROTATE 3\n");
					#pragma acc loop
					for (j=iq+1;j<=n;j++) {			//Sweep columns
						// printf("ip: %d\n",ip);
						// printf("iq: %d\n",iq);
						// printf("j: %d\n",j);
						// printf("a[%d][%d] a[%d][%d]\n",ip,j,iq,j);
						ROTATE(a,n,ip,j,iq,j)
					}
					// printf("ROTATE 4\n");
					#pragma acc loop
					for (j=1;j<=n;j++) {			//Fill eigenvectors
						// printf("ip: %d\n",ip);
						// printf("iq: %d\n",iq);
						// printf("j: %d\n",j);
						// printf("v[%d][%d] v[%d][%d]\n",j,ip,j,iq);
						ROTATE(eigvec,n,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}

		//NO#pragma acc loop	
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];							//Add eigenvalues computed
			eigval[ip]=b[ip];						//Eigenvalues calculated
			z[ip]=0.0;								//Set zero for new iter
		}
		}
	}
	nrerror("Too many iterations in routine jacobi");

}


#undef ROTATE
#undef NRANSI