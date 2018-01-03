subroutine get_byg_ad

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
! Calculate vertical integral of bouyancy gradient (adjoint)           !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  integer(i4)    :: k

      grd%dns(:,:,:) = 0.0

      grd%bx(2:grd%im,:) = grd%bx(2:grd%im,:)/grd%dx(2:grd%im,:)
      grd%by(:,2:grd%jm) = grd%by(:,2:grd%jm)/grd%dy(:,2:grd%jm)

     do k=1,grd%km
      grd%b_x(2:grd%im,:,k) = grd%b_x(2:grd%im,:,k) + grd%bx(2:grd%im,:)*grd%dz(k) * grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k)
      grd%b_y(:,2:grd%jm,k) = grd%b_y(:,2:grd%jm,k) + grd%by(:,2:grd%jm)*grd%dz(k) * grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k)
     enddo

     do k=grd%km,2,-1
      grd%dns(:,2:grd%jm  ,k) = grd%dns(:,2:grd%jm  ,k) +               &
                                grd%b_y(:,2:grd%jm,k)*grd%dz(k)*grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k)
      grd%dns(:,1:grd%jm-1,k) = grd%dns(:,1:grd%jm-1,k) -               &
                                grd%b_y(:,2:grd%jm,k)*grd%dz(k)*grd%msk(:,2:grd%jm,k)*grd%msk(:,1:grd%jm-1,k)
      grd%b_y(:,2:grd%jm,k-1) = grd%b_y(:,2:grd%jm,k-1) + grd%b_y(:,2:grd%jm,k)
      grd%b_y(:,2:grd%jm,k)   = 0.0
     enddo
      grd%dns(:,2:grd%jm  ,1) = grd%dns(:,2:grd%jm  ,1) +               &
                                grd%b_y(:,2:grd%jm,1)*grd%dz(1)*grd%msk(:,2:grd%jm,1)*grd%msk(:,1:grd%jm-1,1)
      grd%dns(:,1:grd%jm-1,1) = grd%dns(:,1:grd%jm-1,1) -               &
                                grd%b_y(:,2:grd%jm,1)*grd%dz(1)*grd%msk(:,2:grd%jm,1)*grd%msk(:,1:grd%jm-1,1)

     do k=grd%km,2,-1
      grd%dns(2:grd%im  ,:,k) = grd%dns(2:grd%im  ,:,k) +               &
                                grd%b_x(2:grd%im,:,k)*grd%dz(k)*grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k)
      grd%dns(1:grd%im-1,:,k) = grd%dns(1:grd%im-1,:,k) -               &
                                grd%b_x(2:grd%im,:,k)*grd%dz(k)*grd%msk(2:grd%im,:,k)*grd%msk(1:grd%im-1,:,k)
      grd%b_x(2:grd%im,:,k-1) = grd%b_x(2:grd%im,:,k-1) + grd%b_x(2:grd%im,:,k)
      grd%b_x(2:grd%im,:,k)   = 0.0
     enddo
      grd%dns(2:grd%im  ,:,1) = grd%dns(2:grd%im  ,:,1) +               &
                                grd%b_x(2:grd%im,:,1)*grd%dz(k)*grd%msk(2:grd%im,:,1)*grd%msk(1:grd%im-1,:,1)
      grd%dns(1:grd%im-1,:,1) = grd%dns(1:grd%im-1,:,1) -               &
                                grd%b_x(2:grd%im,:,1)*grd%dz(k)*grd%msk(2:grd%im,:,1)*grd%msk(1:grd%im-1,:,1)


      grd%tem_ad(:,:,:) = grd%tem_ad(:,:,:) - 0.24*grd%dns(:,:,:) * 9.81/1025. * grd%msr(:,:,:)
      grd%sal_ad(:,:,:) = grd%sal_ad(:,:,:) + 0.74*grd%dns(:,:,:) * 9.81/1025. * grd%msr(:,:,:)

end subroutine get_byg_ad
