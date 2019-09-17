#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

float machineEpsilon(){
    // taking a floating type variable
    float EPS = 0.5;
    float prev_epsilon;
 
    // run until condition satisfy
    while ((1+EPS) != 1)
    {
        // copying value of epsilon into previous epsilon
        prev_epsilon = EPS;
 
        // dividing epsilon by 2
        EPS /=2;
    }
    printf("EPS %f \n", EPS);
    return EPS;
}
void init_mat(int rows, int cols, float **mat){
	int r,c;
	mat[0][0]= 2;
	mat[0][1]= -1;
	mat[0][2]= 0;
	mat[1][0]= -1;
	mat[1][1]= 2;
	mat[1][2]= -1;
	mat[2][0]= 0;
	mat[2][1]= -1;
	mat[2][2]= 2;
}
void print_mat(int rows, int cols, float **mat){
	int r,c;
	printf("\n");
	for (r = 0; r < rows; r++) {
      	for (c = 0; c < cols; c++) 
			printf("%f   ",mat[r][c]);
		printf("\n");
	}
	printf("\n");
}

void max_elem(int *piv_elem,int rows, int cols,float **mat){
	int r,c;
	int max_i = 0;
	int max_j = 1;
	for (r = 0; r < rows-1; r++) 
      	for (c = r+1; c < cols; c++)
      		if(fabs(mat[max_i][max_j])<fabs(mat[r][c])){
      			max_i = r;
      			max_j = c;
    		}
    piv_elem[0] = max_i;
    piv_elem[1] = max_j;
}

float cal_tan(int i,int j,float **mat){
	float num;
	float den;
	float a1;
	float a2;
	float a3;

	num = 2 * (mat[i][j]);
	if(mat[i][i] < mat[i][i])
		num = -num;

	a1 = mat[i][i] - mat[j][j]; 
	a2 = a1*a1;
	a3 = 4 * mat[i][j]*mat[i][j];
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
float cal_sin(float tang, float cose){
	float sino;
	sino = cose*tang;
	return sino;
}

float mat_mult(int rows, int cols,float **mat, float **T, float **mat_temp){
	int i,j,k;

	for (i = 0 ; i < rows ; i++ ){
    	for (j = 0 ; j < cols ; j++ ){
      		mat_temp[i][j] = 0;
        	for (k = 0 ; k < cols ; k++ ){
	            mat_temp[i][j] += T[k][i] * mat[k][j];
			}
		}
	}

	for (i = 0 ; i < rows ; i++ ){
    	for (j = 0 ; j < cols ; j++ ){
      		mat[i][j] = 0;
        	for (k = 0 ; k < cols ; k++ ){
	            mat[i][j] += mat_temp[i][k] * T[k][j];
			}
		}
	}
}

void new_mat(int max_i, int max_j,int rows, int cols,float **mat,float **T,float **mat_temp){
	float tang, cose, sino;
	int c,r;

	tang = cal_tan(max_i,max_j,mat);
	cose = cal_cos(tang);
	sino = cal_sin(cose,tang);

	for (r = 0; r < rows; r++) 
      	for (c = 0; c < cols; c++) 
      		if(r!=c)
      			T[r][c] = 0;
      		else
      			T[r][c] = 1;	

    T[max_i][max_i] = cose;
    T[max_j][max_j] = cose;
    T[max_i][max_j] = -sino; //Termino a eliminar	
    T[max_j][max_i] = sino;		
    

	mat_mult(rows,cols,mat,T,mat_temp);
}

int main(int argc, char const *argv[]){
	int rows = 3, cols = 3;
	int i,j;
	int *piv_elem;
	float **mat, **T, **mat_temp; 
	float sino,cose,tang;
	bool min = false;
	float EPS = machineEpsilon();
	piv_elem=( int *) malloc (2 * sizeof (int));

	mat=( float **) malloc (rows* sizeof (float *));
	for(i=0; i<cols; i++)
		 mat[i] = (float *) malloc(cols * sizeof(float));
	T=( float **) malloc (rows* sizeof (float *));
	for(i=0; i<cols; i++)
		 T[i] = (float *) malloc(cols * sizeof(float));

	mat_temp=( float **) malloc (rows* sizeof (float *));
	for(i=0; i<cols; i++)
		 mat_temp[i] = (float *) malloc(cols * sizeof(float));

	
	init_mat(rows,cols,mat);
	print_mat(rows,cols,mat);

	for (i = 0; min==false; i++){
		printf("iter %d\n",i);
		max_elem(piv_elem,rows,cols,mat);
		new_mat(piv_elem[0],piv_elem[1],rows,cols,mat,T,mat_temp);
		print_mat(rows,cols,mat);

		max_elem(piv_elem,rows,cols,mat);
		if(fabs(mat[piv_elem[0]][piv_elem[1]]) < EPS)
			min=true;
		
	}


	for(i = 0; i < rows; i++) 
		free(T[i]);
	free(T);
	for(i = 0; i < rows; i++) 
		free(mat_temp[i]);
	free(mat_temp);
	for(i = 0; i < rows; i++) 
		free(mat[i]);
	free(mat);
	free(piv_elem);
	return 0;
}