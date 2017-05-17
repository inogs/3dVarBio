
IF(Fortran_COMPILER_NAME STREQUAL "bgxlf90_r")
  foreach(preproc "-Dopt_netcdf" "-Detime=etime_") 
   set(preprocs "${preprocs} -WF,${preproc}")
  endforeach(preproc)
ELSE()
  set(preprocs "-Dopt_netcdf -Detime=etime_")
ENDIF()


# SOURCES --START--
set(KNDSTR  set_knd.f90)

set(OBJSTR  drv_str.f90 cns_str.f90 obs_str.f90 grd_str.f90 bmd_str.f90 eof_str.f90 ctl_str.f90)

set(PHYSOBS bar_mod.f90 bar_mod_ad.f90 div_dmp.f90 div_dmp_ad.f90 get_byg.f90 get_byg_ad.f90 
          get_obs_arg.f90 get_obs_gld.f90 get_obs_xbt.f90 get_obs_sla.f90 
          get_obs_vdr.f90 get_obs_gvl.f90 get_obs_tra.f90 get_obs_trd.f90 
          obs_arg.f90 obs_xbt.f90 obs_gld.f90 obs_vdr.f90 obs_gvl.f90 obs_tra.f90 obs_trd.f90  obs_sla.f90
          obsop_ad.f90 obs_sla_ad.f90 obs_arg_ad.f90 obs_xbt_ad.f90 obs_gld_ad.f90 obs_vdr_ad.f90 obs_gvl_ad.f90 obs_tra_ad.f90
          mod_trj_tl.f90 mod_trj_ad.f90  get_vel_ad.f90 get_vel.f90  
          invrt.f90 invrt_ad.f90 )
set(OBJS routines.f def_nml.f90 def_grd.f90 sav_itr.f90 ini_itr.f90 ini_bmd.f90 rdgrds.f90 rdeofs.f90 
         netcdf_err.f90 get_obs.f90get_obs_sat.f90 
         int_par.f90 obs_vec.f90 def_cov.f90 ini_cfn.f90
         ini_nrm.f90 min_cfn.f90 costf.f90 cnv_ctv.f90 ver_hor.f90 rcfl_x.f90 rcfl_y.f90 veof.f90 
         obsop.f90 obs_sat.f90 resid.f90 res_inc.f90 rdrcorr.f90 mean_rdr.f90     
         obs_trd_ad.f90 obs_sat_ad.f90 veof_ad.f90 ver_hor_ad.f90 rcfl_x_ad.f90 rcfl_y_ad.f90 
         cnv_ctv_ad.f90 cnv_inn.f90 wrt_dia.f90 filename_mod.f90 oceanvar.f90 main.f90)
# SOURCES --END--

add_executable(var_3d ${KNDSTR} ${OBJSTR} ${OBJS})
IF(Fortran_COMPILER_NAME STREQUAL "ifort")
  SET(FLAGS_SPEC "-fltconsistency")
ENDIF()
# LDFLAGS += -L$(LIBDIR) -lf_exit
target_link_libraries(var_3d ${NetCDF_C_LIBRARY} ${FEXIT_LIB})
set_target_properties(var_3d PROPERTIES COMPILE_FLAGS " ${cpp_flags} ${preprocs} ${FLAGS_SPEC}")
include_directories(${NetCDF_INCLUDE_DIR})

install(TARGETS var_3d DESTINATION 3D_VAR)



IF(Fortran_COMPILER_NAME STREQUAL "bgxlf90_r")
 set_SOURCE_FILES_PROPERTIES(routines.f PROPERTIES COMPILE_FLAGS " ${cpp_flags_f} ${preprocs} ${FLAGS_SPEC}")
endif()




