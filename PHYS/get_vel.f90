subroutine get_vel

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
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
! Calculate horizontal velocity from geostrophic formula               !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Bug correction 21.04.2009  thanks to Andrea Storto                   !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  REAL(r8), DIMENSION (grd%im,grd%jm)  :: ud, vd

  integer(i4)    :: k

     do k=grd%km,2,-1
      grd%b_x(2:grd%im,:,k) = ( grd%b_x(2:grd%im,:,k) + grd%b_x(2:grd%im,:,k-1) ) * 0.5
      grd%b_y(:,2:grd%jm,k) = ( grd%b_y(:,2:grd%jm,k) + grd%b_y(:,2:grd%jm,k-1) ) * 0.5
     enddo
      grd%b_x(2:grd%im,:,1) = grd%b_x(2:grd%im,:,1) * 0.5
      grd%b_y(:,2:grd%jm,1) = grd%b_y(:,2:grd%jm,1) * 0.5

     grd%uvl(:,:,:) = 0.0
     grd%vvl(:,:,:) = 0.0

     do k=1,grd%km

      ud(:,:) = 0.0
      vd(:,:) = 0.0

      vd(2:grd%im,:) =   ( (grd%eta(2:grd%im,:)-grd%eta(1:grd%im-1,:))*9.81 + grd%b_x(2:grd%im,:,k) )               &
                           / grd%dx(2:grd%im,:) * grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k) / grd%f(2:grd%im,:)
      ud(:,2:grd%jm) = - ( (grd%eta(:,2:grd%jm)-grd%eta(:,1:grd%jm-1))*9.81 + grd%b_y(:,2:grd%jm,k) )               &
                           / grd%dy(:,2:grd%jm) * grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k) / grd%f(:,2:grd%jm)
      grd%uvl(2:grd%im,1:grd%jm-1,k) =                     &
                         ( ud(2:grd%im,1:grd%jm-1)+ud(1:grd%im-1,1:grd%jm-1)+ud(2:grd%im,2:grd%jm)+ud(1:grd%im-1,2:grd%jm) )*0.25  &
                         * grd%msk(2:grd%im,1:grd%jm-1,k)*grd%msk(1:grd%im-1,1:grd%jm-1,k)
      grd%vvl(1:grd%im-1,2:grd%jm,k) =                     &
                         ( vd(1:grd%im-1,2:grd%jm)+vd(1:grd%im-1,1:grd%jm-1)+vd(2:grd%im,2:grd%jm)+vd(2:grd%im,1:grd%jm-1) )*0.25  &
                         * grd%msk(1:grd%im-1,2:grd%jm,k)*grd%msk(1:grd%im-1,1:grd%jm-1,k)

     enddo


end subroutine get_vel
