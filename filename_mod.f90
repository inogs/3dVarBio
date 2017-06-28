MODULE FILENAMES
implicit none

PUBLIC
character (LEN=1024) :: EOF_FILE_CHL != 'eofs_chl.nc'
character (LEN=1024) :: EOF_FILE_N3N != 'eofs_n3n.nc'
character (LEN=1024) :: EOF_FILE_O2O != 'eofs_o2o.nc'
character (LEN=1024) :: MISFIT_FILE  != 'chl_mis.nc'
character (LEN=1024) :: CORR_FILE    != 'corr.nc'
character (LEN=1024) :: EIV_FILE     != 'eiv.nc'
character (LEN=1024) :: OBS_FILE     != 'obs_1.dat'
character (LEN=1024) :: GRID_FILE    != 'grid1.nc'
!laura
character (LEN=1024) :: RCORR_FILE  != 'chl_rad_corr.nc'
character (LEN=1024) :: ARGO_FILE   != 'argo_mis.dat'
character (LEN=1024) :: ANIS_FILE   != 'gradsal.nc'


CONTAINS
!
!
SUBROUTINE SETFILENAMES
!   !VAR_FILE='DA_static_data/MISFIT/VAR2D/var2D.05.nc'
!
EOF_FILE_CHL = 'eofs_chl.nc'
EOF_FILE_N3N = 'eofs_n3n.nc'
EOF_FILE_O2O = 'eofs_o2o.nc'
MISFIT_FILE  = 'chl_mis.nc'
  CORR_FILE  = 'corr.nc'
   EIV_FILE  = 'eiv.nc'
   OBS_FILE  = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'
  GRID_FILE  = 'grid1.nc'! 'grid'//cgrd//'.nc'
!laura
 RCORR_FILE  = 'chl_rad_corr.nc'
  ARGO_FILE  = 'arg_mis.dat'
  ANIS_FILE  = 'gradsal.nc'

END SUBROUTINE SETFILENAMES






END MODULE FILENAMES
