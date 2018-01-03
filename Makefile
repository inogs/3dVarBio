SHELL = /bin/sh
############################################################################
#
#        Makefile for OceanVar
#
############################################################################
#
#    Copyright 2018 Anna Teruzzi, OGS, Trieste
#
#    This file is part of 3DVarBio.
#    3DVarBio is based on OceanVar (Dobricic, 2006)
#
#    3DVarBio is  free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    3DVarBio is  distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

include compiler.inc

ifndef NETCDF_INC
        export NETCDF_INC=/usr/include
endif

ifndef NETCDF_LIB
        export NETCDF_LIB=/usr/lib
endif
$(info $$NETCDF_INC  = ${NETCDF_INC})
$(info $$NETCDF_LIB  = ${NETCDF_LIB})
$(info $$LIBNCMEDLEV = ${LIBNCMEDLEV})

EXEC = var_3d
LIB  = libvar_3d.a

RM      = rm -f
MV      = mv -f

LIBDEP = libnc-medlevel.a

KNDSTR  =  \
	set_knd.o
OBJSTR  =  \
    filename_mod.o\
		da_params.o\
	drv_str.o\
	cns_str.o\
	obs_str.o\
	grd_str.o\
	eof_str.o\
	ctl_str.o\
	rcfl_mod.o\
	bio_str.o\
	mpi_str.o

PHYSOBS  =  \
	get_obs_sla.o\
	get_obs_arg.o\
	get_obs_xbt.o\
	get_obs_gld.o\
	get_obs_vdr.o\
	get_obs_gvl.o\
	get_obs_tra.o\
	get_obs_trd.o\
	obs_sla.o\
	obs_arg.o\
	obs_xbt.o\
	obs_gld.o\
	obs_vdr.o\
	obs_gvl.o\
	obs_tra.o\
	obs_trd.o\
	obs_sla_ad.o\
	obs_arg_ad.o\
	obs_xbt_ad.o\
	obs_gld_ad.o\
	obs_vdr_ad.o\
	obs_gvl_ad.o\
	obs_tra_ad.o\
	obs_trd_ad.o\
	bar_mod.o\
	get_vel.o\
	div_dmp.o\
	get_byg.o\
	bar_mod_ad.o\
	get_vel_ad.o\
	div_dmp_ad.o\
	get_byg_ad.o\
	mod_trj_ad.o\
	mod_trj_tl.o\
	invrt.o\
	invrt_ad.o

OBJS    =  \
	def_nml.o\
	def_grd.o\
	sav_itr.o\
	rdeofs.o\
    rdrcorr.o\
    mean_rdr.o\
	netcdf_err.o\
	get_obs.o\
	get_obs_arg.o\
	get_obs_sat.o\
	int_par.o\
	obs_vec.o\
	ini_cfn.o\
	cnv_ctv.o\
	ver_hor.o\
	rcfl_x.o\
	rcfl_y.o\
	veof.o\
	obsop.o\
	obs_arg.o\
	resid.o\
	res_inc.o\
	obsop_ad.o\
	obs_arg_ad.o\
	veof_ad.o\
	ver_hor_ad.o\
	rcfl_x_ad.o\
	rcfl_y_ad.o\
	rcfl_y_ad_init.o\
	rcfl_y_init.o\
	cnv_ctv_ad.o\
	cnv_inn.o\
	wrt_dia.o\
	clean_mem.o\
	mpi_utils.o\
	readAnisotropy.o\
	bio_mod.o\
	bio_mod_ad.o\
	readBioStat.o\
	wrt_bio_stat.o\
	costf.o\
	obs_sat.o\
	bio_conv.o\
	bio_conv_ad.o\
	obs_sat_ad.o\
	readGrid.o\
	def_cov.o\
	tao_minimizer.o\
    oceanvar.o

MAINEXE = main.o

.SUFFIXES: .f90

all:  $(EXEC) $(LIB)
	@echo $(EXEC) is compiled

install:	$(EXEC)
		cp -p $(EXEC) $(INSTDIR)

$(EXEC) : $(LIBDEP)	$(KNDSTR) $(OBJSTR) $(OBJS) $(MAINEXE)
	$(LD)  -o $(EXEC) $(OBJSTR) $(OBJS) $(MAINEXE) $(LDFLAGS)

$(LIB)  :       $(KNDSTR) $(OBJSTR) $(OBJS)
	ar -r $(LIB) $(KNDSTR) $(OBJSTR) $(OBJS)

tao_minimizer.o: tao_minimizer.f90
	$(CPP) -I$(PETSC_INC) $*.f90 > cpp.$*.f90 ; $(F90) -I$(PETSC_INC) $(FFLAGS) cpp.$*.f90  ; $(MV) cpp.$*.o $*.o

mpi_utils.o: mpi_utils.f90
	$(CPP) -I$(PETSC_INC) $*.f90 > cpp.$*.f90 ; $(F90) -I$(PETSC_INC) $(FFLAGS) cpp.$*.f90  ; $(MV) cpp.$*.o $*.o

.DEFAULTS:
.f90.o :
	$(CPP) $*.f90 > cpp.$*.f90 ; $(F90) $(FFLAGS) cpp.$*.f90  ; $(MV) cpp.$*.o $*.o

.f.o :
	$(CPP) $*.f > cpp.$*.f ; $(F77) $(FFLAGS) cpp.$*.f  ; $(MV) cpp.$*.o $*.o

libf_exit.a :
	cd $(LIBFEXIT) && $(MAKE)

libnc-medlevel.a :
		cd $(LIBNCMEDLEV) && $(MAKE)

clean:
	$(RM) *.o *.mod cpp.* *.L
	cd $(LIBNCMEDLEV) && $(MAKE) erase
	cd ..

erase:
	$(RM) *.o *.mod cpp.* *.L $(EXEC)
	cd $(LIBNCMEDLEV) && $(MAKE) erase
	cd ..
