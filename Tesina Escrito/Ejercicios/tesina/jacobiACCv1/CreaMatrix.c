#include <stdio.h>
#define NP 128
void main()
{
	float m[NP][NP];
	int i,j,k;
	for ( i = 0; i < NP; ++i)
	{
		for (j = i,k=NP/2; j < NP; ++j,k--)
		{
			m[i][j]=m[j][i]= k;
		}
	}

	for ( i = 0; i < NP; ++i){
		//printf("\n");
	for ( j = 0; j < NP; ++j)
	{
		printf("%.1f,", m[i][j]);
	}

	printf("\n");
}
}