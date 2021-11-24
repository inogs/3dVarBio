subroutine res_inc

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
! Initialise for adjoint calculations                                  !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use grd_str
 use obs_str
 use bio_str

 implicit none

 if (drv%multiv .eq. 0) then
  if (drv%chl_assim .eq. 1) then
    grd%chl_ad(:,:,:) = 0.0 ! OMP
  end if
  
  if (drv%nut .eq. 1) then
    if (bio%n3n .eq. 1) &
      grd%n3n_ad(:,:,:) = 0.0
    if (bio%o2o .eq. 1) &
      grd%o2o_ad(:,:,:) = 0.0
  endif

 else if(drv%multiv .eq.1) then
    grd%chl_ad(:,:,:) = 0.0 ! OMP
    grd%n3n_ad(:,:,:) = 0.0
 endif
 
 obs%gra(:) = obs%amo(:) / obs%err(:) ! OMP

end subroutine res_inc
