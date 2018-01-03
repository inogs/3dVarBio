  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2018 Anna Teruzzi, OGS, Trieste                         !
  !                                                                          !
  !    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                          !
  !                                                                          !
  !    3DVarBio is  free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    3DVarBio is  distributed in the hope that it will be useful,           !
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !---------------------------------------------------------------------------

MODULE FILENAMES
implicit none

PUBLIC
character (LEN=1024) :: EOF_FILE    != 'eofs.nc'
character (LEN=1024) :: MISFIT_FILE != 'chl_mis.nc'
character (LEN=1024) :: CORR_FILE   != 'corr.nc'
character (LEN=1024) :: EIV_FILE    != 'eiv.nc'
character (LEN=1024) :: OBS_FILE    != 'obs_1.dat'
character (LEN=1024) :: GRID_FILE   != 'grid1.nc'
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
   EOF_FILE = 'eofs.nc'
MISFIT_FILE = 'chl_mis.nc'
  CORR_FILE = 'corr.nc'
   EIV_FILE = 'eiv.nc'
   OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'
  GRID_FILE = 'grid1.nc'! 'grid'//cgrd//'.nc'
!laura
 RCORR_FILE = 'chl_rad_corr.nc'
  ARGO_FILE = 'arg_mis.dat'
  ANIS_FILE = 'gradsal.nc'

END SUBROUTINE SETFILENAMES






END MODULE FILENAMES
