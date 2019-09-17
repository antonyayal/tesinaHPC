void eigsrt(float d[], float *v, int n)
{
	int k,j,i;
	float p;
	#pragma acc loop
	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j*n+i];
				v[j*n+i]=v[j*n+k];
				v[j*n+k]=p;
			}
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
