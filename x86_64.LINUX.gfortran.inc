FC      = gfortran $(OPENMP_FLAG)
CC      = gcc
F90     = gfortran $(OPENMP_FLAG)
F77     = gfortran
CPP     = gfortran -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
MPFC    = mpif90
MPCC    = mpicc
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _
OPT_FLTCONSISTENCY =

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

# NETCDF_INC = $(NETCDFF_INC)
# NETCDF_LIB = $(NETCDFF_LIB)
FFLAGS  = -O2 -I$(NETCDF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -lnetcdff
AR      = ar

# Debug options##########
# FFLAGS = -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -fbacktrace -c $(OPT_FLTCONSISTENCY)
#########################

