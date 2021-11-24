subroutine obsop_ad
  
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
  ! Apply observational operators - adjoint
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use obs_str
  use drv_str
  use mpi_str
  
  implicit none
  
  INTEGER(i4)  :: ierr
  
  obs%k = 0

  ! ---
  ! ARGO observations
  if (drv%argo_obs .eq. 1) &
    call obs_arg_ad

  ! ---
  ! Observations of chlorophyll
  if(drv%sat_obs .eq. 1) &  
    call obs_sat_ad

  ! ---
  ! Apply biological repartition of the chlorophyll
  if((drv%chl_assim .eq. 1) .or. (drv%multiv .eq. 1)) &
    call bio_conv_ad

  call MPI_Barrier(Var3DCommunicator, ierr)
  
end subroutine obsop_ad
