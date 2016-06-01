FC      = ifort $(OPENMP_FLAG)
CC      = icc
F90     = ifort $(OPENMP_FLAG)
F77     = ifort
CPP     = ifort -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
MPFC    = mpiifort
MPCC    = mpiicc
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _
OPT_FLTCONSISTENCY =

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

FFLAGS  = -g -O2 -I$(NETCDF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -g -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -lnetcdff -lnetcdf
AR      = ar

# Debug options##########
# FFLAGS =  -O0 -g -traceback -fp-stack-check -check bounds -fpe0 -c $(OPT_FLTCONSISTENCY)
#########################

