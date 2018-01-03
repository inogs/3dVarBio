subroutine get_byg

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
! Calculate vertical integral of bouyancy gradient                     !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  INTEGER(i4)    :: k


! ---
! Assume that the density is a simple linear function of temperature and salinity

     grd%dns(:,:,:) = (-0.24*grd%tem(:,:,:) + 0.74*grd%sal(:,:,:))*9.81/1025. * grd%msk(:,:,:)

! ---
! Bouyancy force
      grd%b_x(:,:,:) = 0.0
      grd%b_x(2:grd%im,:,1) =                                               &
             (grd%dns(2:grd%im,:,1)-grd%dns(1:grd%im-1,:,1))*grd%dz(1)*grd%msk(2:grd%im,:,1)*grd%msk(1:grd%im-1,:,1)
      grd%b_y(:,:,:) = 0.0
      grd%b_y(:,2:grd%jm,1) =                                               &
             (grd%dns(:,2:grd%jm,1)-grd%dns(:,1:grd%jm-1,1))*grd%dz(1)*grd%msk(:,2:grd%jm,1)*grd%msk(:,1:grd%jm-1,1)
     do k=2,grd%km
      grd%b_x(2:grd%im,:,k) = grd%b_x(2:grd%im,:,k-1) +                     &
                         (grd%dns(2:grd%im,:,k)-grd%dns(1:grd%im-1,:,k))*grd%dz(k)*grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k)
      grd%b_y(:,2:grd%jm,k) = grd%b_y(:,2:grd%jm,k-1) +                     &
                         (grd%dns(:,2:grd%jm,k)-grd%dns(:,1:grd%jm-1,k))*grd%dz(k)*grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k)
     enddo

! ---
! Verical integral of bouyancy force
      grd%bx(:,:) = 0.0
      grd%by(:,:) = 0.0
     do k=1,grd%km
      grd%bx(2:grd%im,:) = grd%bx(2:grd%im,:) + grd%b_x(2:grd%im,:,k)*grd%dz(k) * grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k)
      grd%by(:,2:grd%jm) = grd%by(:,2:grd%jm) + grd%b_y(:,2:grd%jm,k)*grd%dz(k) * grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k)
     enddo

      grd%bx(2:grd%im,:) = grd%bx(2:grd%im,:)/grd%dx(2:grd%im,:)
      grd%by(:,2:grd%jm) = grd%by(:,2:grd%jm)/grd%dy(:,2:grd%jm)

end subroutine get_byg
