#ifndef LANCZOS_H
#define LANCZOS_H

#include<stdbool.h>
#ifndef SEQ_PARA
#define SEQ_PARA
#endif

#ifdef GSL_BLAS_FONCTION
#undef SEQ_PARA
#endif

#define ForD(i,debut,fin,step) for(i = debut; i < fin; i +=step)

 /* Algorithme de Lanczos */
#ifdef SEQ_PARA
extern void Lanczos(const double** A,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m);

extern void LanczosSparse(const int *row,const int * col, const double * val,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m);

#endif

#ifdef GSL_BLAS_FONCTION
extern void Lanczos(const double * const A,double ** Vm, double** T, double * r0,uint* itermax, uint nrowA, uint ncolA, uint nrowVm, uint m);
#endif

 /* Test si vm est une base */

extern bool IsBase(const double ** const Vm, uint m, uint nrowVm);

/* Somme de deux vecteur */
extern double * SomVect(const double *  u, const double * v, double * som, uint taille);

 /* Difference de deux vecteur */
extern double * DiffVect(const double * u, const double * v, double * diff,uint taille);

 /* Difference mult de deux vecteurs */
extern double * FMDIFF(double * u,const double * v,const double alpha,uint taille);


 /* Multiplication par un scalaire */
double * MultVectScal(const double scalaire, const double * u,double * scla_u, uint taille);

double * MultVectScalRef(const double scalaire,double * u,uint taille);
 /* Produit scalaire de deux vecteurs */
extern double ProdScal(const double * u, const double * v, uint taille);

 /* Produit Matrice vecteur */
extern double * ProdMatVect(const double * * A, const double * X, double * Y,uint rowA, uint colA, uint sizeX, uint sizeY);

/* Produit Matrice vecteur diff */
extern double * ProdMatVectDiff(const double ** A, const double * X,const double alpha,const double * Y,double * R,uint nrowA, uint ncolA, uint sizeX, uint sizeY,uint sizeR);

/* Produit Sparse Matrice vecteur */
extern double * ProdMatVectDiffSparse(const int * row,const int * col,const double * val, const double * X,const double alpha,const double * Y,double * R,uint nrowA, uint ncolA, uint sizeX, uint sizeY,uint sizeR);

 /* Produit Matrice Matrice */
extern double ** ProdMatMat(const double ** A, const double ** B, double ** C,uint rowA,uint colA,uint rowB,uint colB,uint rowC, uint colC);



#endif
