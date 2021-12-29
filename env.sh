
#!/bin/bash

# Les chemins vers les bibliotheques
GCC_DIR_9_2_0=/gpfs/softs/spack/opt/spack/linux-centos7-haswell/gcc-4.8.5/gcc-9.2.0-k4debouieljbsurbrr2w755davwl6xsu

export GCC_INCLUDE_DIR_9_2_0=${GCC_DIR_9_2_0}/lib/gcc/x86_64-pc-linux-gnu/9.2.0/include
export GCC_LIB_DIR_9_2_0=${GCC_DIR_9_2_0}/lib

OPENMPI_DIR=/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/gcc-9.2.0/openmpi-4.0.2-uhldlns6vmbx54du6p2h5fyimtrye6dx
export MPI_INCLUDE_DIR_4_0_2=${OPENMP_DIR}/include
export MPI_LIB_DIR_4_0_2=${OPENMP_DIR}/lib

GSL_DIR=/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/gcc-9.2.0/gsl-2.6-qxs6txqmfqb2qccdo6pm5uwz65zhmzpc
export GSL_INCLUDE_DIR_2_6=${GSL_DIR}/include/gsl
export GSL_LIB_DIR_2_6=${GSL_DIR}/lib

OPENBLAS_DIR=/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/gcc-9.2.0/openblas-0.3.8-br7jxs72pbbccb6xsbdcqwgks3pdfvcx
export BLAS_INCLUDE_DIR_0_3_8=${OPENBLAS_DIR}/include
export BLAS_LIB_DIR_0_3_8=${OPENBLAS_DIR}/lib

# Chargement des modules
export MODULES_LIST="openmpi/4.0.2/gcc-9.2.0 gsl/2.6/gcc-9.2.0 openblas/0.3.8/gcc-9.2.0 gcc/9.2.0/gcc-4.8.5"
module load ${MODULES_LIST}

#num_erreur=$?

#echo "La derniere commande "
#echo $num_erreur

#if ["$num_errer" -ne "0"]
#then
# exit $num_erreur
#fi

# Affichage des informations importantes
echo ""
echo "**************Les differents modules chargés voir dans MODULES_LIST*******"

for m in ${MODULES_LIST} 
do
    echo $m
done


echo ""
echo "************La version du compilateur************* "
gcc --version
echo "which gcc"
which gcc
