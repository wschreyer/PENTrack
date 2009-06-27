#include <malloc.h>
#include <stdio.h>
#include "nrutil.h"
#include "main.h"

void nrerror(char *error_text){	
	//void exit(int c);

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(-1);
}


/*
long double *vector(int nl,int nh)
//int nl,nh;
{
	long double *v;

	v=(long double *)malloc((unsigned) (nh-nl+1)*sizeof(long double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}
*/
int *ivector(int nl, int nh)
//int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

long double *dvector(int nl, int nh)
//int nl,nh;
{
	long double *v;

	v=(long double *)malloc((unsigned) (nh-nl+1)*sizeof(long double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

long double ****viertensor(int n1,int nh1,int n2,int nh2,int n3,int nh3,int n4,int nh4)
{
	int i=0,k=0;
	long double ****tensor;
	
	tensor = (long double ****) malloc((unsigned) (nh1-n1+1)*sizeof(long double ***));
	if(!tensor) nrerror("allocation failure 1 in viertensor()");
		
	tensor -= n1;
	
	for(i=n1;i<=nh1;i++){
		tensor[i] = (long double ***) malloc((unsigned) (nh2-n2+1)*sizeof(long double **));
		if(!tensor[i]) nrerror("allocation failure 2 in viertensor()");
		tensor[i] -= n2;
		for(k=n2;k<=nh2;k++)
			tensor[i][k] = matrix(n3,nh3,n4,nh4);
	}
	
	return tensor;
	
}

void free_viertensor(long double ****tensor,int n1,int nh1,int n2,int nh2,int n3,int nh3,int n4,int nh4)
{
	int i=0,k=0;
	for(i=n1;i<=nh1;i++){
		for(k=n2;k<=nh2;k++)
			free_matrix(tensor[i][k],n3,nh3,n4,nh4);
		free(tensor[i]+n2);
	}
	free(tensor+n1);
}

long double **matrix(int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
	int i;
	long double **m;

	m=(long double **) malloc((unsigned) (nrh-nrl+1)*sizeof(long double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(long double *) malloc((unsigned) (nch-ncl+1)*sizeof(long double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}


long double **save_matrix(int nrl, int nrh, int ncl, int nch)
{
	int i=0,j=0;
	long double **m;

	m=(long double **) malloc((unsigned) (nrh-nrl+1)*sizeof(long double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(long double *) malloc((unsigned) (nch-ncl+1)*sizeof(long double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
		for(j=1;j<=nch;j++)
			m[i][j]=0.0;
	}
	return m;
}

long double **dmatrix(int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
	int i;
	long double **m;

	m=(long double **) malloc((unsigned) (nrh-nrl+1)*sizeof(long double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(long double *) malloc((unsigned) (nch-ncl+1)*sizeof(long double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

long double **submatrix(long double **a, int oldrl, int oldrh, int oldcl,int oldch,
      int newrl, int newcl)
/*long double **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;*/
{
	int i,j;
	long double **m;

	m=(long double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(long double*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(long double *v, int nl, int nh)

{
	free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)

{
	free((char*) (v+nl));
}

void free_dvector(long double *v, int nl, int nh)

{
	free((char*) (v+nl));
}

void free_matrix(long double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(long double **m, int nr1, int nrh, int nc1, int nch)
{
	int i;
	for(i=nrh;i>=nr1;i--) free(m[i]+nc1);
	free(m+nr1);
}

void free_imatrix(int **m, int nr1, int nrh, int nc1, int nch)
{
	int i;
	for(i=nrh;i>=nr1;i--) free(m[i]+nc1);
	free(m+nr1);
}

void free_submatrix(long double **b, int nrl,int nrh, int ncl,int nch)
{
	free((char*) (b+nrl));
}

long double **convert_matrix(long double *a, int nrl, int nrh, int ncl, int nch)
{
	int i,j,nrow,ncol;
	long double **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (long double **) malloc((unsigned) (nrow)*sizeof(long double*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}

void free_convert_matrix(long double **b,int nrl, int nrh,int ncl,int nch)

{
	free((char*) (b+nrl));
}
