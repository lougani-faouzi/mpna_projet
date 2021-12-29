#ifndef MATRICE_HPP
#define MATRICE_HPP

#define For(i,debut,fin,step) for(i = debut; i < fin; i +=step)

double** A1_matrice();
double** A3_matrice(int nrow,int ncol,int debut, int fin);
double** A9_matrice(int nrow,int ncol,int debut, int fin);
double** AM_n_matrice(int nrow, int ncol);
double** MatTest();


int AM_n_matrice_sparse(int *row,int *col,double *val,int nrow,int ncol);
int A9_matrice_sparse(int *row,int *col,double *val,int nrow,int ncol);


double* A3_matrice_blas(int nrow,int ncol,int debut, int fin);
double* A9_matrice_blas(int nrow,int ncol,int debut, int fin);
void VecteurAleatoire(double * x,int size);

void PrintMatrice(double* A,int nrow, int ncol);
#endif
