# Choix du compilateur
CC =gcc
DEBUG=-g

#macros
#MACROS_BLAS=-DHAVE_INLINE -DBIB_GSL_BLAS
#MATRICE_SPARSE=-DSPARSEMATRIX
#VERSION_PARALLEL=-DPARALLEL
#VERSION_BLAS=-DGSL_BLAS_FONCTION
# Flags
FLAGS_OPTIMISATION=-g -fopt-info-vec -march=native $(MACROS_BLAS) $(VERSION_BLAS) $(VERSION_PARALLEL)
CLIBS=-fopenmp -lm -latomic -lpthread
Warnings= -W -Wall -pedantic

# Les repertoires des codes
ACUTUAL_DIR :=.
SRCDIR=$(ACUTUAL_DIR)/src
INCLDIR=$(ACUTUAL_DIR)/include
OBJDIR=$(ACUTUAL_DIR)/obj
BINDIR=$(ACUTUAL_DIR)/bin

# Recuperation des fichiers : *.c, *.h
srcs= $(wildcard $(SRCDIR)/*.c)
incls := $(subst $(SRCDIR),$(INCLDIR),$(subst .c,.h,$(srcs)))
objs :=$(subst $(SRCDIR),$(OBJDIR),$(subst .c,.o,$(srcs)))

#----------EXECUTABLE---------
EXEC=$(BINDIR)/main

# MISE A JOUR DES VARIABLES D'ENVIRONNEMENT

##-----Pour blas
INCLUDE_BIB_BLAS=-I$(BLAS_INCLUDE_DIR_0_3_8)
LIB_BIB_BLAS=-L$(BLAS_LIB_DIR_0_3_8) -lopenblas

##---Pour gsl (blas avec matrice creuse)
INCLUDE_BIB_GSL=-I$(GSL_INCLUDE_DIR_2_6)
LIB_BIB_GSL=-L$(GSL_LIB_DIR_2_6) -lgsl -lgslcblas

##---Pour mpi
#INCLUDE_BIB_MPI=$(shell mpicc --showme:compile)
#LIB_BIB_MPI=$(shell mpicc --showme:link)

##---Pour GCC
INCLUDE_LIB_GCC=
LIB_BIB_GCC=$(CLIBS)


INCLUDE_BIB=$(INCLUDE_LIB_GCC) $(INCLUDE_BIB_GSL) $(INCLUDE_BIB_BLAS) $(INCLUDE_BIB_MPI)
LIB_BIB=$(LIB_BIB_GSL) $(LIB_BIB_BLAS) $(LIB_BIB_MPI) $(LIB_BIB_GCC)

# PREMIERE REGLE
.PHONY:all
all: makedir $(EXEC) 

#------COMPILATION DU CODE LANCZOS --------

$(EXEC):main.o $(objs)
	$(CC) $(DEBUG) $(Warnings) $(FLAGS_OPTIMISATION) -I$(INCLDIR) $(INCLUDE_BIB) $^ -o $@ $(LIB_BIB)

main.o : main.c $(incls)
	$(CC) $(DEBUG) $(Warnings) $(FLAGS_OPTIMISATION) -I$(INCLDIR) $(INCLUDE_BIB) -c $< -o $@ $(LIB_BIB)

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(INCLDIR)/%.h
	$(CC) $(DEBUG) $(Warnings) $(FLAGS_OPTIMISATION) -I$(INCLDIR) $(INCLUDE_BIB) -c $< -o $@ $(LIB_BIB)
	echo $(MODULES_LIST)



# Netoyage et/ou creation des repertoires
.PHONY:clean makedir ans
clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(BINDIR)/$(EXEC)
	rm -f main.o

makedir:
	mkdir -p $(SRCDIR)
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(INCLDIR)

ans:$(EXEC)
	$(EXEC)
