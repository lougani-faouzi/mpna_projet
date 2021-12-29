#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>
#include<stdbool.h>
#include <math.h>

/* ----Commande du proc -----------*/
#define SEQ_PARA
#define DENSEMATRIX

#ifdef SPARSEMATRIX
#undef DENSEMATRIX
#endif

#ifdef GSL_BLAS_FONCTION
#undef SEQ_PARA
#endif
/* --------Inclusion des header internes------------*/
#include"include/lanczos.h"
#include"include/matrice.h"

/*--------Iclusion des header externes--------------*/
#ifdef BIB_GSL_BLAS
#include<gsl_cblas.h>
#include<gsl_spmatrix.h>
#include<gsl_vector.h>
#include<gsl_spblas.h>
#include <gsl/gsl_eigen.h>
#endif

#ifdef PARALLEL
#include "omp.h"
#define nb_thread 4
#endif
/*-------------Partie code --------------------------*/
#define gettime(t) clock_gettime(CLOCK_MONOTONIC_RAW, t)
#define get_sub_seconde(t) (1e-9*(double)t.tv_nsec)
/** return time in second
*/
double get_elapsedtime(void)
{
  struct timespec st;
  int err = gettime(&st);
  if (err !=0) return 0;
  return (double)st.tv_sec + get_sub_seconde(st);
}


#define SIZE 10000
#define M 5

void Print(double ** A, int nrow,int ncol,int debut,int fin);
void Print2(double* A,int nrow, int ncol,int debut,int fin);

#ifdef BIB_GSL_BLAS
void ValeursPropres(const double** const A,const double** const H,const size_t nrow, const size_t ncol, const size_t size);
#endif



int main(int argc, char * argv[])
{  
  int nrow = SIZE; int ncol = SIZE;
  int m = M;
  
  if( 2 <= argc) { nrow = atoi(argv[1]); ncol = nrow ;}
  if( 3 <= argc) {m = atoi(argv[2]);}  
#ifdef PARALLEL
  int nbt = nb_thread;
  if(4 <= argc) {nbt = atoi(argv[3]);}
  omp_set_num_threads(nbt);
#endif
  int debut = 0; int fin = ncol-1;
  
#ifdef GSL_BLAS_FONCTION
  double *A = A9_matrice_blas(nrow,ncol,debut,fin);
#endif

  /* Cas sequentiel ou parallel */
#ifdef SEQ_PARA
/* matrice dense */
  #ifdef DENSEMATRIX
 
   double ** A = A9_matrice(nrow,ncol,debut,fin); 
  #endif
/* Sparse matrix */
  #ifdef SPARSEMATRIX

  int * RowA=NULL; int * ColA=NULL; double * ValA=NULL;
  
  size_t bytes = ceil((double)nrow* 1./64.)*64*3*sizeof(int);
  int test = posix_memalign((void**)&RowA,64,bytes);
  assert(!test);
  
  size_t bytes_64 = ceil((double)nrow* 1./64.)*64*3*sizeof(int);
  int test1 = posix_memalign((void**)&ColA,64,bytes_64);
  assert(!test1);

   size_t bytes_64_1 = ceil((double)nrow* 1./64.)*64*3*sizeof(double);
   int test3 = posix_memalign((void**)&ValA,64,bytes_64_1);
  assert(!test3);
  
  uint nnz = A9_matrice_sparse(RowA,ColA,ValA,nrow,ncol);
  printf("NNZ = %d \n",nnz);
  #endif
#endif
  
  double ** Vm = (double**) malloc(m*sizeof(double*));
  int j = 0;
  ForD(j,0,m,1) Vm[j] = (double*) malloc(nrow*sizeof(double));  
  double ** T =  (double**) malloc(m*sizeof(double*));
  ForD(j,0,m,1) T[j] = malloc(nrow*sizeof(double));
  
  /*-------- Initialisation des vecteurs--------- */  
  double * r0= (double*) malloc(nrow*sizeof(double));
  VecteurAleatoire(r0,nrow);  
  
  uint MaxIter = 0;
  
  /*----------Appel de l'algorithme---------------*/
  double t0 = 0., t1 = 0., duration = 0.;
  t0 = get_elapsedtime();
  /* Cas matrice Dense */
#ifdef DENSEMATRIX
  Lanczos(A,Vm,T,r0,&MaxIter,nrow,ncol, nrow,m);
#endif
  /* Cas matrice creuse */
#ifdef SPARSEMATRIX
  LanczosSparse(RowA,ColA,ValA,Vm,T,r0,&MaxIter,nrow,ncol, nrow,m);
#endif
  t1 = get_elapsedtime();
   
  printf("Matrice A\n");
  /* Print2(A,nrow,ncol,debut,fin); */
  int i = 0;
  /* printf("Matrice Vm\n"); */
  /* Print(Vm,m,ncol,0,m-1); */

  printf("Matrice Tm %d,\n",MaxIter);
  /* Print(T,MaxIter,MaxIter,0,MaxIter-1); */

  /* Calcul des valeurs propres */
  #ifdef BIB_GSL_BLAS
  /* ValeursPropres(A,T,nrow,ncol,MaxIter); */
  #endif

  // Pretty print
  duration = (t1 - t0);
  fprintf(stdout, "Performance results: \n");
  fprintf(stdout, "  Time: %lf s\n", duration);
  
  
  /*------------Liberation de memoire-------------*/
#ifdef GSL_BLAS_FONCTION
  free(A);
#endif

#ifdef SEQ_PARA

#ifdef DENSEMATRIX
ForD(i,0,m,1) free(A[i]);
free(A);
#endif
 
#ifdef SPARSEMATRIX
 free(RowA); free(ColA); free(ValA);
#endif

#endif
 
  ForD(i,0,m,1) free(Vm[i]);
  free(Vm);
  ForD(i,0,m,1) free(T[i]);
  free(T);
  free(r0);

  return 0;
}



void Print2(double* A,int nrow, int ncol,int debut,int fin)
{
  int i = 0, j=0;
  int size = fin - debut +1;
  assert(fin < nrow);
  ForD(i,0,size,1)
    {
      ForD(j,0,ncol,1)
	{
	  printf("%lf  ",A[i*ncol+j]);
	}
      printf("\n");
    }
}

void Print(double** A,int nrow, int ncol,int debut,int fin)
{
  int i = 0, j=0;
  int size = fin - debut +1;
  assert(fin < nrow);
  ForD(i,0,size,1)
    {
      ForD(j,0,ncol,1)
	{
	  printf("%lf  ",A[i][j]);
	}
      printf("\n");
    }
}

#ifdef BIB_GSL_BLAS
void ValeursPropres(const double** const A,const double** const H,const size_t nrow, const size_t ncol, const size_t size)
{
  gsl_matrix * Agsl = gsl_matrix_alloc(nrow,ncol);
  /* gsl_matrix_memcpy(Hm,H); */
  size_t i = 0;
  size_t j = 0;
  /* Initialisation des tableaux */
  ForD(i,0,nrow,1)
  {
    ForD(j,0,ncol,1)
      {
	double alpha = A[i][j];
	gsl_matrix_set(Agsl,i,j,alpha);
      }
  }
  assert(nrow==ncol);
  /* Pour la deuxieme matrice */
  gsl_matrix * Hm = gsl_matrix_alloc(size,size);
  /* Pour i = 0 */
  double alpha = H[0][0]; double beta = H[0][1]; double beta_inf = beta;
  /* double norm = sqrt(alpha *alpha + beta * beta); */
  /* alpha = alpha / norm; beta = beta/norm; */
  gsl_matrix_set(Hm,0,0,alpha);gsl_matrix_set(Hm,0,1,beta);
  /* Pour i = size-1 */
  alpha = H[size-1][size-1]; beta = H[size-1][size-2];  
  gsl_matrix_set(Hm,size-1,size-1,alpha); gsl_matrix_set(Hm,size-1,size-2,beta);
  /* Pour le reste de i */
  ForD(i,1,size-1,1)
  {
    alpha = H[i][i];
    beta = H[i][i+1];
    /* norm = sqrt(alpha *alpha + 2*beta * beta); */
    gsl_matrix_set(Hm,i,i,alpha);
    gsl_matrix_set(Hm,i,i+1,beta);
    gsl_matrix_set(Hm,i+1,i,beta_inf);
    beta_inf = beta;
  }
  /* La valeur et vecteur propre pour Hm */
  gsl_vector * ValPropreHm = gsl_vector_alloc(size);
  gsl_matrix * VectPropreHm = gsl_matrix_alloc(size,size);

  /* Les valeurs et vecteurs propre pour A */
  gsl_vector * ValPropreA = gsl_vector_alloc(nrow);
  gsl_matrix * VectPropreA = gsl_matrix_alloc(nrow,ncol);

  /* Pour les paramettres */
  gsl_eigen_symmv_workspace * wHm = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv_workspace * wA = gsl_eigen_symmv_alloc(nrow);
  
  gsl_eigen_symmv(Hm,ValPropreHm,VectPropreHm,wHm);
  gsl_eigen_symmv(Agsl,ValPropreA,VectPropreA,wA);

   /* On libere de la memoire */
  gsl_eigen_symmv_free(wHm);
  gsl_eigen_symmv_free(wA);
  
  /* Pour la matrice A */
  gsl_eigen_symmv_sort(ValPropreA,VectPropreA,GSL_EIGEN_SORT_VAL_DESC);
  {
    uint i;
    ForD(i,0,nrow,1)
    {
      double ValPropreA_i= gsl_vector_get(ValPropreA,i);
      gsl_vector_view VectPropreA_i = gsl_matrix_column(VectPropreA,i);
      printf("Valeur propre Pour A = %g ", ValPropreA_i);
      /* printf("Vecteur propre pour A = \n"); */
      /* gsl_vector_fprintf(stdout,&VectPropreA_i.vector, "%g"); */
    }
    printf("\n");
  }
  
  /* Pour la deuxieme matrice */
  gsl_eigen_symmv_sort(ValPropreHm,VectPropreHm,GSL_EIGEN_SORT_VAL_DESC);
   {
    uint i;
    ForD(i,0u,size,1)
    {
      double ValPropreHm_i = gsl_vector_get(ValPropreHm,i);
      gsl_vector_view VectPropreHm_i = gsl_matrix_column(VectPropreHm,i);
      printf("Valeur propre Pour Hm = %g\n", ValPropreHm_i);
      /* printf("Vecteur propre pour Hm = \n"); */
      /* gsl_vector_fprintf(stdout,&VectPropreHm_i.vector, "%g"); */
    }
  }

   /* On libere de la memeoire */
   gsl_vector_free(ValPropreHm);
   gsl_vector_free(ValPropreA);
   
   gsl_matrix_free(VectPropreHm);
   gsl_matrix_free(VectPropreA);
  
}
#endif
