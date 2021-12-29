#include <stdlib.h>
#include <assert.h>
#include "../include/matrice.h"
#include <math.h>
#include <stdio.h>


double** A1_matrice()
{
  double a =11111111,b=9090909,c=10891089,d=8910891;
  double e =11108889,f=9089091,g=10888911,h=8909109;
  double** A = (double**) malloc(8*sizeof(double*));
  uint i = 0;
  For(i,0,8u,1) A[i] = (double*) malloc(8*sizeof(double*));

  A[0][0] = a; A[0][1]=-b;A[0][2]=-c;A[0][3]=d; A[0][4]=-e;A[0][5]=f; A[0][6]=g; A[0][7]=-h;
  A[1][0] =-b; A[1][1]=a;A[1][2]=d;A[1][3]=-c; A[1][4]=f;A[1][5]=-e; A[1][6]=-h; A[1][7]=g;
  A[2][0] =-c; A[2][1]=d;A[2][2]=a;A[2][3]=-b; A[2][4]=g;A[2][5]=-h; A[2][6]=-e; A[2][7]=f;
  A[3][0] =d; A[3][1]=-c;A[3][2]=-b;A[3][3]=a; A[3][4]=-h;A[3][5]=g; A[3][6]=f; A[3][7]=-e;
  A[4][0] =-e; A[4][1]=f;A[4][2]=g;A[4][3]=-h; A[4][4]=a;A[4][5]=-b; A[4][6]=-c; A[4][7]=d;
  A[5][0] =f; A[5][1]=-e;A[5][2]=-h;A[5][3]=g; A[5][4]=-b;A[5][5]=a; A[5][6]=d; A[5][7]=-c;
  A[6][0] =g; A[6][1]=-h;A[6][2]=-e;A[6][3]=f; A[6][4]=-c;A[6][5]=d; A[6][6]=a; A[6][7]=-b;
  A[7][0] =-h; A[7][1]=g;A[7][2]=f;A[7][3]=-e; A[7][4]=d;A[7][5]=-c; A[7][6]=-b; A[7][7]=a;

  return A;
}
/* Matrice A3 */
double** A3_matrice(int nrow, int ncol,int debut,int fin)
{
  assert(nrow==ncol);
  assert(fin < nrow);
  assert(debut <= fin);
  int size = fin-debut+1;
  size_t bytes = size*sizeof(double*);
  double ** A = (double**) malloc(bytes);
  int i = 0;
  For(i,0,size,1)
    {
      A[i] = (double*) malloc(ncol*sizeof(double));
    }
 
 /* Initialisation */
 
 for(int i = 0; i < size; i++)
   {
     int I = i+debut;
     for(int j = 0; j < I; j++)
     {
       A[i][j] = nrow+1-I;               
     }

     for(int j = I; j < ncol; j++)
     {
       A[i][j] =nrow+1 - j;               
     }
   }
  

 return A;
}

/* Matrice A9 */

double** A9_matrice(int nrow, int ncol,int debut, int fin)
{
  assert(nrow==ncol);
  assert(fin < nrow);
  assert(debut <= fin);
  int size = fin-debut+1;
  size_t bytes = (size)*sizeof(double*);
  double ** A = (double**) malloc(bytes);
  int i = 0;
  For(i,0,size,1)
    {
      A[i] = (double*) malloc(ncol*sizeof(double));
    }

  double a = 11111111, b = 9090909;

  for (int i = 0; i <size; i++)
    {
      if (debut+i != 0) A[i][i+debut-1] = b;
      A[i][debut+i] = a;
      if (i+debut !=nrow-1)  A[i][i+debut+1] = b;
    }
  
  return A;
}

/* SPARSE A9 */
int A9_matrice_sparse(int *RowA,int *ColA,double *ValA,int nrow,int ncol)
{
  assert(nrow==ncol);

  double a = 11111111, b = 9090909;
  
  int i = 0; 
  RowA[0] = 0; RowA[1] = 2;
  ColA[0] = 0; ColA[1] = 1;
  ValA[0] = a; ValA[1] = b;
  uint nnz=2;
  
  For(i,1,nrow-1,1)
  {
    RowA[i+1] = RowA[i] + 3;
    ColA[nnz] = i-1; ColA[nnz+1] = i; ColA[nnz +2] = i+1;
    ValA[nnz] = b; ValA[nnz+1] = a; ValA[nnz +2] = b;
    nnz += 3; 
  }
  RowA[nrow] = RowA[nrow-1]+2;
  ColA[nnz] = nrow-2; ColA[nnz+1] = nrow-1;
  ValA[nnz] = b;      ValA[nnz+1] = a;

  return nnz;
}

/* Matrice An */

double** AM_n_matrice(int nrow, int ncol)
{
  size_t bytes = nrow*sizeof(double*);
  double ** A = (double**) malloc(bytes);
  int i = 0;
  For(i,0,nrow,1)
    {
      A[i] = (double*) malloc(ncol*sizeof(double));
    }

  A[0][0] = 1; A[0][1] = -1e-1;
  A[nrow-1][ncol-2] = 1e-1; A[nrow-1][ncol-1] = nrow;

  for(int i = 1; i < nrow-1; i++)
    {
      A[i][i] = i+1;
      A[i][i-1] = 1e-1;
      A[i][i+1] = -1e-1;
    }
  return A;
}

double ** MatTest()
{
  int nrow = 203;
  size_t bytes = nrow*sizeof(double*);
  double ** A = (double**) malloc(bytes);

  int i = 0;
  For(i,0,nrow,1)
    {
      A[i] = (double*) malloc(nrow*sizeof(double));
    }
  
  double s = 0.0;

  for(int i = 0; i < nrow-2; i++)
    {
      A[i][i] = s;
      s += 0.01;
    }
  A[201][201] = 2.5; A[202][202] = 3.;
  return A;
}

/* Partie blas */
double* A3_matrice_blas(int nrow,int ncol,int debut, int fin)
{
  assert(nrow==ncol);
  assert(fin < nrow);
  assert(debut <= fin);
  int offset = fin-debut+1;
  double * A = NULL;

  size_t bytes_64 = ceil((double)ncol* 1./64.)*64*offset*sizeof(double);
  uint test = posix_memalign((void**)&A,64,bytes_64);
  assert(!test);
 
 /* Initialisation */
 
 for(int i = 0; i < offset; i++)
   {
     int I = i+debut;
     assert(I<ncol);
     for(int j = 0; j < I; j++)
     {
       A[i*ncol +j] = nrow+1-I;               
     }

     for(int j = I; j < ncol; j++)
     {
       A[i*ncol+j] = nrow+1 - j;               
     }
   }

  return A;
}



  
double* A9_matrice_blas(int nrow,int ncol,int debut, int fin)
{
  assert(nrow==ncol);
  assert(fin < nrow);
  assert(debut <= fin);
  int offset = fin-debut+1;
  double * A = NULL;

  size_t bytes_64 = ceil((double)ncol* 1./64.)*64*offset*sizeof(double);
  uint test = posix_memalign((void**)&A,64,bytes_64);
  assert(!test);

  double a = 11111111, b = 9090909;

  for (int i = 0; i <offset; i++)
    {
      if (debut+i != 0) A[i*ncol+i+debut-1] = b;
      A[i*ncol+debut+i] = a;
      if (i+debut !=nrow-1)  A[i*ncol+i+debut+1] = b;
    }
  
  return A;
}


void VecteurAleatoire(double * x,int size)
{
  srand(42);
  int i = 0;
  For(i,0,size,1)
    {
      x[i] = (double)rand()/(double)RAND_MAX;
    }
}


void PrintMatrice(double* A,int nrow, int ncol)
{
  int i = 0, j=0;
  
  For(i,0,nrow,1)
    {
      For(j,0,ncol,1)
	{
	  printf("%g  ",A[i*ncol+j]);
	}
      printf("\n");
    }
}

