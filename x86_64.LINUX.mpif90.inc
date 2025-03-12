FC      = mpif90 $(OPENMP_FLAG)
CC      = mpicc
F90     = mpif90 $(OPENMP_FLAG)
F77     = mpif77
CPP     = mpif90 -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
MPFC    = mpif90
MPCC    = mpicc
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

PETSC_LD_FLAGS := $(if $(PETSC_LIB), -L$(PETSC_LIB),)

# NETCDF_INC = $(NETCDFF_INC)
# NETCDF_LIB = $(NETCDFF_LIB)

PNETCDF_INC=/m100_scratch/userexternal/gbolzon0/pnetcdf/include/
PNETCDF_LIB=/m100_scratch/userexternal/gbolzon0/pnetcdf/lib

# FFLAGS  = -ffree-line-length-none -O2 -I$(NETCDF_INC) $(OPT_FLTCONSISTENCY) -c
FFLAGS  = -ffree-line-length-none -O2 -I$(NETCDF_INC) -I$(NETCDFF_INC) -I$(PNETCDF_INC) -c
CFLAGS  = -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -L$(NETCDFF_LIB) $(PETSC_LD_FLAGS) -L$(PNETCDF_LIB) -lnetcdff -lnetcdf -lpetsc -lpnetcdf
AR      = ar

# Debug options##########
# FFLAGS = -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -fbacktrace -c $(OPT_FLTCONSISTENCY)
#########################

