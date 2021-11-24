subroutine obsop
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
  !                                                                          !
  !    This file is part of OceanVar.                                          !
  !                                                                          !
  !    OceanVar is free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    OceanVar is distributed in the hope that it will be useful,           !
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !---------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply observational operators   
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------

  
  use set_knd
  use obs_str
  use drv_str
  use mpi_str
  
  implicit none
  
  INTEGER(i4) :: ierr

  call MPI_Barrier(Var3DCommunicator, ierr)

  ! ---
  ! Apply biological repartition of the chlorophyll
  if((drv%chl_assim .eq. 1) .or. (drv%multiv .eq. 1)) &
    call bio_conv

  ! ---
  ! Observations by ARGO floats
  if (drv%argo_obs .eq. 1) &
    call obs_arg
  
  ! ---
  ! Observations of chlorophyll
  if(drv%sat_obs .eq. 1) &
    call obs_sat

  call MPI_Barrier(Var3DCommunicator, ierr)
  
end subroutine obsop
