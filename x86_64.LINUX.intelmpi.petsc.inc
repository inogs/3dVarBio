FC      = mpiifort $(OPENMP_FLAG)
CC      = mpiicc
F90     = mpiifort $(OPENMP_FLAG)
F77     = mpiifort
CPP     = mpiifort -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
MPFC    = mpiifort
MPCC    = mpiicc
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _
OPT_FLTCONSISTENCY =

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

PREPROC = -D_USE_MPI
CPP += $(PREPROC)

# NETCDF_INC = $(NETCDFF_INC)
# NETCDF_LIB = $(NETCDFF_LIB)

# FFLAGS  = -ffree-line-length-none -O2 -I$(NETCDF_INC) -I$(NETCDFF_INC) $(OPT_FLTCONSISTENCY) -c
FFLAGS  = -g -O2 -I$(NETCDF_INC) -I$(NETCDFF_INC) -I$(PNETCDF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -L$(NETCDFF_LIB) -L$(PETSC_LIB) -L$(PNETCDF_LIB) -lnetcdff -lnetcdf -lpetsc -lpnetcdf
AR      = ar

# Debug options##########
# FFLAGS = -O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -fbacktrace -c $(OPT_FLTCONSISTENCY)
#########################

