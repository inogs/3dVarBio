subroutine int_par

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
! Calculate interpolation parameters                                   !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use obs_str

 implicit none

#ifdef __FISICA
! ----
! Load SLA observations
  call int_par_sla

! ----
! Load ARGO observations
  call int_par_arg

! ----
! Load XBT observations
  call int_par_xbt

! ----
! Load glider observations
  call int_par_gld

! ----
! Load observations of Argo trajectories
  call int_par_tra

! ----
! Load observations of trajectories of surface drifters
  call int_par_trd

! ----
! Load observations of drifter velocities
  call int_par_vdr

! ----
! Load observations of glider velocities
  call int_par_gvl

#endif
! ----
! Load observations of chlorophyll
  call int_par_chl

end subroutine int_par
