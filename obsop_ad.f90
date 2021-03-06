subroutine obsop_ad
  
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
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply observational operators - adjoint
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
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
  if(drv%bio_assim .eq. 1) &
    call bio_conv_ad

  call MPI_Barrier(Var3DCommunicator, ierr)
  
end subroutine obsop_ad
