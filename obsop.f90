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

 implicit none

#ifdef __FISICA
! ---
! Satellite observations of SLA
  call obs_sla

! ---
! Observations by ARGO floats
  call obs_arg

! ---
! Observations by XBT profiles
  call obs_xbt

! ---
! Observations by gliders
  call obs_gld

! ---
! Observations of Argo trajectories
  if(tra%no.gt.0) call obs_tra

! ---
! Observations of trajectories of surface drifters
  if(trd%no.gt.0) call obs_trd

! ---
! Observations of velocities by drifters
  call obs_vdr

! ---
! Observations of velocities by gliders
  call obs_gvl

#endif
! ---
! Observations of chlorophyll
  call obs_chl

end subroutine obsop
