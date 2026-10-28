# Compiler settings for the build inside an OPA chain.
# Compilers (Intel oneAPI ifx/icx with Intel MPI) and libraries (NetCDF,
# PnetCDF, PETSc) are taken from the spack environment of the chain, in
# $(OPA_HOME)/HOST/$(OPA_HOSTNAME)/spack_env; no module has to be loaded.

ifeq ($(strip $(OPA_HOME)),)
$(error OPA_HOME is not defined: this file can be used only inside an OPA chain)
endif
ifeq ($(strip $(OPA_HOSTNAME)),)
$(error OPA_HOSTNAME is not defined: this file can be used only inside an OPA chain)
endif

SPACK_VIEW = $(OPA_HOME)/HOST/$(OPA_HOSTNAME)/spack_env/.spack-env/view

# The Intel MPI wrappers find their installation through I_MPI_ROOT and the
# compilers (ifx, icx) through PATH
export I_MPI_ROOT := $(SPACK_VIEW)/mpi/latest
export PATH := $(SPACK_VIEW)/oneapi-compilers/compiler/latest/bin:$(I_MPI_ROOT)/bin:$(PATH)

NETCDF_INC  = $(SPACK_VIEW)/include
NETCDF_LIB  = $(SPACK_VIEW)/lib
NETCDFF_INC = $(SPACK_VIEW)/include
NETCDFF_LIB = $(SPACK_VIEW)/lib
PNETCDF_INC = $(SPACK_VIEW)/include
PNETCDF_LIB = $(SPACK_VIEW)/lib
PETSC_INC   = $(SPACK_VIEW)/include
PETSC_LIB   = $(SPACK_VIEW)/lib

FC      = mpiifx $(OPENMP_FLAG)
CC      = mpiicx
F90     = mpiifx $(OPENMP_FLAG)
F77     = mpiifx
CPP     = mpiifx -cpp -E
# -static-intel links the Intel runtime libraries statically: libimf.so needs
# libintlc.so.5, that the RUNPATH of the executable does not cover, so in a
# clean environment it would be found only if another library (e.g. PETSc)
# happened to load it
LD      = $(FC) $(OPENMP_FLAG) -static-intel
MPFC    = mpiifx
MPCC    = mpiicx
MPLD    = $(MPFC)

FORTRAN_UNDERSCORE = _
OPT_FLTCONSISTENCY =

LIBFEXIT  = ./libfexit/
LIBNCMEDLEV =./libnc-medlevel

PREPROC =
CPP += $(PREPROC)

FFLAGS  = -g -O2 -I$(NETCDF_INC) -I$(NETCDFF_INC) -I$(PNETCDF_INC) $(OPT_FLTCONSISTENCY) -c
CFLAGS  = -O2 -I$(NETCDF_INC)

LDFLAGS += -L$(NETCDF_LIB) -L$(NETCDFF_LIB) -L$(PETSC_LIB) -L$(PNETCDF_LIB) -lnetcdff -lnetcdf -lpetsc -lpnetcdf
# The libraries of the spack view are found at run time through the paths
# written in the executable, so it runs also in a clean environment.
# The Intel MPI wrappers write them as RUNPATH (--enable-new-dtags, that can
# not be overridden), so LD_LIBRARY_PATH, if set, takes precedence over them.
LDFLAGS += -Wl,-rpath,$(SPACK_VIEW)/lib
AR      = ar
