/*
 * creaMatrix.c
 *
 *
 * José Antonio Ayala Barbosa
 * Oct, 2018
 *
 * 
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>


#define NP 16
int main(int argc, char const *argv[]){

	float m[NP][NP];
	int i,j,k;
	char nombreArchivo[100];
	FILE *archivo;
	char tempChar[100];

	for ( i = 0; i < NP; ++i)
	{
		for (j = i,k=NP/2; j < NP; ++j,k--)
		{
			m[i][j]=m[j][i]= k;
		}
	}

	strcpy(nombreArchivo, "matrix");
    sprintf(tempChar, "%d", NP);
    strcat(nombreArchivo, tempChar);
    strcat(nombreArchivo, ".txt");

	if (fopen(nombreArchivo, "r") == NULL){ // Si no existe archivo crea un archivo con los numeros
			printf("Creating File\n");
			archivo = fopen(nombreArchivo, "w"); // abre archivo de escritura
			fprintf(archivo,"%d",NP);
			fprintf(archivo,"%c",'\n');
		   	
		   	for (i = 0; i < NP; i++){
   				for (j = 0; j < NP; j++){
   					
		    		fprintf(archivo,"%.1f ",m[i][j]); //Escribe numeros aleatorios 
		    	}
		    	fprintf(archivo,"%c",'\n');
			}
		    fclose(archivo); // Cierra el archivo
		
			printf("File created\n");
	}

	printf("\n");
	return 0;
}