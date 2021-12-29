# cisseAly_khadirSelma_louganiFaouzi_codeSourceLanczos

* Ce code s'execute uniquement sur le cluster RUCHE et pas en local ,une fois on a reserver un noeud:
* il faut charger l'environement avec le script suivant:
```
source env.sh

```
* Pour compiler on tappe la commande :
```
make

```
* Pour executer la version sequentielle dense il faut laisser toute la partie macros commentee Dans le makefile
```
#macros
#MACROS_BLAS=-DHAVE_INLINE -DBIB_GSL_BLAS
#MATRICE_SPARSE=-DSPARSEMATRIX
#VERSION_PARALLEL=-DPARALLEL
#VERSION_BLAS=-DGSL_BLAS_FONCTION
```
* Puis Taper (pour 10000=taille de la matrice A et nb valeur propres=5) apartir du dossier principal:

```
./bin/main 10000 5

```
* Pour executer la version parralele dense : decomenter VERSION_PARALLEL=-DPARALLEL
executer la commande avec 4 thread OpenMp

```
./bin/main 10000 5 4

```
* Pour executer la version blas dense: decomenter les 2 lignes MACROS_BLAS=-DHAVE_INLINE -DBIB_GSL_BLAS et VERSION_BLAS=-DGSL_BLAS_FONCTION
apres tappez:

```
./bin/main 10000 5

```
* Pour executer la version sequentielle creuse : decomenter MATRICE_SPARSE=-DSPARSEMATRIX
apres tappez:
```
./bin/main 10000 5

```

* Pour executer la version parralele creuse : decomenter les 2 lignes VERSION_PARALLEL=-DPARALLEL et MATRICE_SPARSE=-DSPARSEMATRIX
executer la commande avec 4 thread OpenMp

```
./bin/main 10000 5 4

```

* Pour supprimer les executables:
```
make clean
```


