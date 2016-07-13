FC      = mpifort $(OPENMP_FLAG)
CC      = mpicc
F90     = mpifort $(OPENMP_FLAG)
F77     = mpifort
CPP     = mpifort -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
MPFC    = mpifort
MPCC    = mpicc
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _
OPT_FLTCONSISTENCY =

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

PREPROC = -D_USE_MPI
CPP += $(PREPROC)

# NETCDF_INC = $(NETCDFF_INC)
# NETCDF_LIB = $(NETCDFF_LIB)

FFLAGS  = -ffree-line-length-none -O2 -I$(NETCDF_INC) $(OPT_FLTCONSISTENCY) -c
# FFLAGS  = -O2 -I$(NETCDF_INC) -I$(NETCDFF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -L$(PETSC_LIB) -lnetcdff -lnetcdf -lpetsc
AR      = ar

# Debug options##########
# FFLAGS = -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -fbacktrace -c $(OPT_FLTCONSISTENCY)
#########################

