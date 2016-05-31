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

FFLAGS  = -O0 -g -traceback -fp-stack-check -check bounds -fpe0 -I$(NETCDF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -O0 -g -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -lnetcdff
AR      = ar

