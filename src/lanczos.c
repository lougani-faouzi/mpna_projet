#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include <math.h>

#include"../include/lanczos.h"



#define SEQ_PARA

/*--------- Partie gsl_blas--------------- */
#ifdef GSL_BLAS_FONCTION
#undef SEQ_PARA
#undef PARALLEL
#include<gsl_cblas.h>
#include<gsl_vector.h>
#include <gsl_matrix.h>
#endif

/* ------Partie OMP----------- */
#ifdef PARALLEL
#include "omp.h"
#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline") //Optimization flags
#endif

/* Somme de deux vecteur */
inline double * SomVect(const double *  u, const double * v, double * som, uint taille)
{
  uint i = 0;
#ifdef PARALLEL
  #pragma omp parallel for
#endif
  ForD(i,0,taille,1)  som[i] = u[i] + v[i];
  return som;
}

/* Difference de deux vecteur */
inline double * DiffVect(const double *  u, const double * v, double * diff,uint taille)
{
  uint i = 0;
#ifdef PARALLEL
#pragma omp parallel for
#endif
  ForD(i,0,taille,1)  diff[i] = u[i] - v[i];
  return diff;
}

 /* Difference, mult vecteurs */
inline double * FMDIFF(double * u,const double * v,const double alpha,uint taille)
{
  uint i = 0;
#ifdef PARALLEL
  #pragma omp parallel for
#endif
  ForD(i,0,taille,1)  u[i] = u[i] - alpha*v[i];
  return u;
}

/* Multiplication par un scalaire */
double * MultVectScal(const double scalaire,const double * u,double * scla_u, uint taille)
{
  uint i = 0;
#ifdef PARALLEL
  #pragma omp parallel for
#endif
  ForD(i,0,taille,1) scla_u[i] = scalaire * u[i];
  return scla_u;
}

double * MultVectScalRef(const double scalaire,double * u,uint taille)
{
  uint i = 0;
#ifdef PARALLEL
#pragma omp parallel for
#endif
  ForD(i,0,taille,1) u[i] = scalaire * u[i];
  return u;
}

/*--------- Produit scalaire de deux vecteurs--------- */
double ProdScal(const double * u, const double * v, uint taille)
{
  double s = 0.0;
  uint i = 0;
#ifdef PARALLEL
  #pragma omp parallel
  #pragma omp for reduction(+:s)
#endif
  ForD(i,0,taille,1) s += u[i]*v[i];       
  return s;
}
 
/* Produit Matrice vecteur */
double * ProdMatVect(const double ** A, const double * X, double * Y,uint rowA, uint colA, uint sizeX, uint sizeY)
{
  uint i = 0;
  uint j = 0;
  // Test sur les taille
  assert(rowA==sizeY);
  assert(colA==sizeX);

#ifdef PARALLEL
#pragma omp parallel
#pragma omp for
#endif
  ForD(i,0,rowA,1)
    {
      double s = 0.;
      ForD(j,0,colA,1)
	{
	  s += A[i][j]*X[j];
	}
      Y[i] = s;
    }

  return Y;
}

/*------- Produit Matrice vecteur diff mult---------- */
inline double * ProdMatVectDiff(const double ** A, const double * X,const double alpha,const double * Y,double * R,uint nrowA, uint ncolA, uint sizeX, uint sizeY,uint sizeR)
{
  
  uint i = 0;
  uint j = 0;

  // Test sur les taille
  assert(nrowA==sizeY);
  assert(ncolA==sizeX);
  assert(sizeR==sizeY);

  if(R==NULL) exit(1);/* R = (double*) malloc(sizeR *sizeof(double)); */

#ifdef PARALLEL
#pragma omp parallel for
#endif
  ForD(i,0,nrowA,1)
  {
      double s = 0.;
      ForD(j,0,ncolA,1)
      {
	  s += A[i][j]*X[j];
      }
      R[i] = s- alpha*Y[i];
    }

  return R;
}


/* Pour matrice creuse */
double * ProdMatVectDiffSparse(const int * row,const int * col,const double * val, const double * X,const double alpha,const double * Y,double * R,uint nrowA, uint ncolA, uint sizeX, uint sizeY,uint sizeR)
{
  int i = 0;
  int j = 0;
  // Test sur les taille
  assert(nrowA==sizeY);
  assert(ncolA==sizeX);
  assert(sizeR==sizeY);
  
#ifdef PARALLEL
#pragma omp parallel for
#endif
  
  ForD(i,0,nrowA,1)
  {
    int debut = row[i];
    int fin = row[i+1];
    double s = 0;
    ForD(j,debut,fin,1)
    {
      uint col_index = col[j];
      s += val[j] * X[col_index];
    }
    R[i] = s - alpha * Y[i];
  }
  return R;
}


/* Produit Matrice Matrice */
double ** ProdMatMat(const double ** A, const double ** B, double ** C,uint rowA,uint colA,uint rowB,uint colB,uint rowC, uint colC)
{
  uint i = 0;
  uint j = 0;
  uint k = 0;
  // Test sur les taille
  assert(rowA==rowC);
  assert(colA==rowB);
  assert(colB==colC);  

  ForD(i,0,rowA,1)
    {
      ForD(k,0,colB,1)
	{
	  double aik = A[i][k];
	  ForD(j,0,colA,1)
	    {
	      C[i][j] += aik * B[k][j];
	    }
	}  
    }
  return C;
}


/*---------- Algorithme de Lanczos---------------*/

#ifdef SEQ_PARA

void Lanczos(const double ** A,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m)
{
  
#ifdef PARALLEL
   printf("---Version Parallele activée dans lanczos-- \n");
#endif

   #ifndef PARALLEL
   printf("---Version Sequentielle dans lanczos-- \n");
   #endif
  assert(nrowVm == ncolA);
  assert(m <= nrowA);
  /* Varibles */
  double beta_j = 0.;
  double beta_jp1 = 0.;
  double alpha_j = 0.;
  double const epselon = 1e-17;
  double norm_r0 = sqrt(ProdScal(r0,r0,nrowVm));
  Vm[0] = MultVectScal(1./norm_r0,r0,Vm[0],nrowVm);
  
  
  uint j = 0u;
  double * w =NULL;
  uint bytes_mult_64 = ceil((double)nrowVm* 1./64.)*64*sizeof(double);
  uint test = posix_memalign((void**)&w,64,bytes_mult_64);
  assert(test==0);

  ForD(j,0,m,1)
    {
      uint i = 0;
      if(j!=0) i = j-1;     
      w = ProdMatVectDiff(A,Vm[j],beta_j,Vm[i],w,nrowA,ncolA,nrowVm,nrowVm,nrowVm);
      alpha_j = ProdScal(w,Vm[j],nrowVm);
      w = FMDIFF(w,Vm[j],alpha_j,nrowVm);
      beta_jp1 = sqrt(ProdScal(w,w,nrowVm));
      beta_j = beta_jp1;
      
      T[j][j] = alpha_j;
      if(j+1 < m) T[j][j+1] = beta_jp1;
      if(j!=0) {T[j][j-1] = T[j-1][j];/* printf("T[j-1][j]=%g",T[j-1][j]); */}

      if(beta_jp1<epselon){*itermax=j+1; break; }
      if(j+1<m) Vm[j+1]= MultVectScal(1./beta_jp1,w,Vm[j+1],nrowVm);
      *itermax= j+1;
      /* if(j+1<m) printf("is ortho? %g, norm= %g \n",ProdScal(Vm[j],Vm[j+1],nrowVm),ProdScal(Vm[j+1],Vm[j+1],nrowVm)); */
    }
  free(w);
}


/*--------- Pour la matrice creuse-----------*/
void LanczosSparse(const int *RowA,const int * ColA, const double * ValA,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m)
{
    
#ifdef PARALLEL
   printf("---Version Sparse Parallele activée dans lanczos-- \n");
#endif

#ifndef PARALLEL
   printf("---Version Sparse Sequentielle dans lanczos-- \n");
#endif
   
  assert(nrowVm == ncolA);
  assert(m <= nrowA);
  /* Varibles */
  double beta_j = 0.;
  double beta_jp1 = 0.;
  double alpha_j = 0.;
  double const epselon = 1e-17;
  double norm_r0 = sqrt(ProdScal(r0,r0,nrowVm));
  Vm[0] = MultVectScal(1./norm_r0,r0,Vm[0],nrowVm);
  
  
  uint j = 0u;
  double * w =NULL;
  size_t bytes_mult_64 = ceil((double)nrowVm* 1./64.)*64*sizeof(double);
  uint test = posix_memalign((void**)&w,64,bytes_mult_64);
  assert(test==0);

  ForD(j,0,m,1)
    {
      uint i = 0;
      if(j!=0) i = j-1;     
      w = ProdMatVectDiffSparse(RowA,ColA,ValA,Vm[j],beta_j,Vm[i],w,nrowA,ncolA,nrowVm,nrowVm,nrowVm);
      alpha_j = ProdScal(w,Vm[j],nrowVm);
      w = FMDIFF(w,Vm[j],alpha_j,nrowVm);
      beta_jp1 = sqrt(ProdScal(w,w,nrowVm));
      beta_j = beta_jp1;
      
      T[j][j] = alpha_j;
      if(j+1 < m) T[j][j+1] = beta_jp1;
      if(j!=0) {T[j][j-1] = T[j-1][j];/* printf("T[j-1][j]=%g",T[j-1][j]); */}

      if(beta_jp1<epselon){*itermax=j+1; break; }
      if(j+1<m) Vm[j+1]= MultVectScal(1./beta_jp1,w,Vm[j+1],nrowVm);
      *itermax= j+1;
      /* if(j+1<m) printf("is ortho? %g, norm= %g \n",ProdScal(Vm[j],Vm[j+1],nrowVm),ProdScal(Vm[j+1],Vm[j+1],nrowVm)); */
    }
  free(w);
}

#endif



 /*---------- Pour le cas blas ---------*/
#ifdef GSL_BLAS_FONCTION
void Lanczos(const double * A,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m)
{
  printf("---Version Blas activé dans lanczos-- \n");
  assert(nrowVm == ncolA);
  assert(m <= nrowA);
  /* Varibles */
  double beta_j = 0.;
  double beta_jp1 = 0.;
  double alpha_j = 0.;
  double const epselon = 1e-17;
  
  double norm_r0 = cblas_dnrm2(nrowVm,r0,1);
  cblas_dcopy(nrowVm,r0,1,Vm[0],1);
  cblas_dscal(nrowVm,1./norm_r0,Vm[0],1);
  
  uint j = 0u;
  double * w =NULL;
  uint bytes_mult_64 = ceil((double)nrowVm* 1./64.)*64*sizeof(double);
  uint test = posix_memalign((void**)&w,64,bytes_mult_64);
  assert(test==0);

  ForD(j,0,m,1)
    {
      uint i = 0;
      if(j!=0) i = j-1;
      
      cblas_dgemv(CblasRowMajor,CblasNoTrans,nrowA,ncolA,1.,A,nrowVm,Vm[j],1,0.,w,1);
      alpha_j = cblas_ddot(nrowVm,w,1,Vm[j],1);
      cblas_daxpy(nrowVm,-1*beta_j,Vm[i],1,w,1);
      cblas_daxpy(nrowVm,-1*alpha_j,Vm[j],1,w,1);
      beta_jp1 = cblas_dnrm2(nrowVm,w,1);

       beta_j = beta_jp1;
      /* Matrice de Ritz */      
      T[j][j] = alpha_j;
      if(j+1 < m) T[j][j+1] = beta_jp1;
      if(j!=0) {T[j][j-1] = T[j-1][j];/* printf("T[j-1][j]=%g",T[j-1][j]); */}
      
      if(beta_jp1<epselon){*itermax=j+1; break; }
      if(j+1<m)
	{
	  cblas_dcopy(nrowVm,w,1,Vm[j+1],1);
	  cblas_dscal(nrowVm,1./beta_jp1,Vm[j+1],1);
	}
      *itermax= j+1;
      /* if(j+1<m) printf("is ortho? %g, norm= %g \n",cblas_ddot(nrowVm,Vm[j],1,Vm[j+1],1),cblas_ddot(nrowVm,Vm[j+1],1,Vm[j+1],1)); */
    }
  free(w);
}
#endif

