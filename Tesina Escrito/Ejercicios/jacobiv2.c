#include <math.h>
#define NRANSI
#include "nrutil.h"

#define ROTATE(a,i,j,k,l) \
		EPS=a[i][j];\
		temp=a[k][l];\
		a[i][j]=EPS-s*(temp+EPS*tau);\
		a[k][l]=temp+s*(EPS-temp*tau);

void jacobil(float **a, int n, int ndm, float **eigvec, float eigenvalues[],  int *nrot){

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

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) { 						//Initializating identity matrix.
		for (iq=1;iq<=n;iq++) eigvec[ip][iq]=0.0;
		eigvec[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {						
		b[ip]=eigenvalues[ip]=a[ip][ip];			//Initializating b and eigvec to the diagonal af a.
		z[ip]=0.0;									//This vector will accumulate terms of the form tapq as in equation. (creo que old vector)
	}
	*nrot=0;									
	for (i=1;i<=50;i++) {
		sm=0.0;					
		for (ip=1;ip<=n-1;ip++) {					//Sum off-diagonal upper triangle elements.
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}

		if (sm == 0.0) {							//Finish jacobi,
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}

		if (i < 4)									//First three steps.
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;							
		
		for (ip=1;ip<=n-1;ip++) {

			for (iq=ip+1;iq<=n;iq++) {
				EPS=100.0*fabs(a[ip][iq]); 			//Epsilon of the machine

				if (i > 4 							//After 4 sweeps, skip rotation if the off-diagonal element is too small. 
					&& (float)(fabs(eigval[ip])+EPS) == (float)fabs(eigval[ip])
					&& (float)(fabs(eigval[iq])+EPS) == (float)fabs(eigval[iq]))
					a[ip][iq]=0.0;					//As element is too small, = 0
				
				else if (fabs(a[ip][iq]) > tresh){ 	//If still can be computed.
					
					temp=eigval[iq]-eigval[ip];

					if ((float)(fabs(temp)+EPS) == (float)fabs(temp))		/******CHECAR******/
						t=(a[ip][iq])/temp;							//tanΘ = 1/(2Θ)
					
					else {
						theta=0.5*temp/(a[ip][iq]);					//Θ = (aqq-app)/2apq
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));	//tanΘ = ±1/(|Θ|+√(Θ^2+1))
						if (theta < 0.0) t = -t;					//plus sign is used if app≥aqq, else minus sign
					}

					c=1.0 / sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					temp=t*a[ip][iq];
					z[ip] -= temp;
					z[iq] += temp;
					eigenvalues[ip] -= temp; 		//Row
					eigenvalues[iq] += temp;		//Column
					a[ip][iq]=0.0;
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}

		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];							//Add eigenvalues computed
			eigenvalues[ip]=b[ip];					//Eigenvalues calculated
			z[ip]=0.0;								//Set zero for new iter
		}
	}
	nrerror("Too many iterations in routine jacobi");
}


#undef ROTATE
#undef NRANSI