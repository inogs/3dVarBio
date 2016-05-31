MODULE FILENAMES
implicit none

PUBLIC
character (LEN=1024) :: EOF_FILE    != 'eofs.nc'
character (LEN=1024) :: MISFIT_FILE != 'chl_mis.nc'
character (LEN=1024) :: CORR_FILE   != 'corr.nc'
character (LEN=1024) :: EIV_FILE    != 'eiv.nc'
character (LEN=1024) :: OBS_FILE    != 'obs_1.dat'
character (LEN=1024) :: GRID_FILE   != 'grid1.nc'




CONTAINS
!
!
SUBROUTINE SETFILENAMES
!   !VAR_FILE='DA_static_data/MISFIT/VAR2D/var2D.05.nc'
!
   EOF_FILE = 'eofs.nc'
MISFIT_FILE = 'chl_mis.nc'
  CORR_FILE = 'corr.nc'
   EIV_FILE = 'eiv.nc'
   OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'
  GRID_FILE = 'grid1.nc'! 'grid'//cgrd//'.nc'

END SUBROUTINE SETFILENAMES






END MODULE FILENAMES
