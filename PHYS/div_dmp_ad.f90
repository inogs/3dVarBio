subroutine div_dmp_ad

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
! Divergence damping (adjoint)
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str

  implicit none

  INTEGER(i4)    :: k, kdiv
  REAL(r8), DIMENSION (grd%im,grd%jm)  :: div


      div(:,:) = 0.0

     do k=grd%km,1,-1


      do kdiv = 1,100

       div(1:grd%im-1,2:grd%jm-1) = div(1:grd%im-1,2:grd%jm-1) + grd%vvl(1:grd%im-1,2:grd%jm-1,k)*0.2 * grd%adxdy**2    &
                                  / grd%dy(1:grd%im-1,2:grd%jm-1)                                                       &
                                  * grd%msk(1:grd%im-1,2:grd%jm-1,k)*grd%msk(1:grd%im-1,1:grd%jm-2,k)
       div(1:grd%im-1,1:grd%jm-2) = div(1:grd%im-1,1:grd%jm-2) - grd%vvl(1:grd%im-1,2:grd%jm-1,k)*0.2 * grd%adxdy**2    &
                                  / grd%dy(1:grd%im-1,2:grd%jm-1)                                                       &
                                  * grd%msk(1:grd%im-1,2:grd%jm-1,k)*grd%msk(1:grd%im-1,1:grd%jm-2,k)
       div(2:grd%im-1,1:grd%jm-1) = div(2:grd%im-1,1:grd%jm-1) + grd%uvl(2:grd%im-1,1:grd%jm-1,k)*0.2 * grd%adxdy**2    &
                                  / grd%dx(2:grd%im-1,1:grd%jm-1)                                                       &
                                  * grd%msk(2:grd%im-1,1:grd%jm-1,k)*grd%msk(1:grd%im-2,1:grd%jm-1,k)
       div(1:grd%im-2,1:grd%jm-1) = div(1:grd%im-2,1:grd%jm-1) - grd%uvl(2:grd%im-1,1:grd%jm-1,k)*0.2 * grd%adxdy**2    &
                                  / grd%dx(2:grd%im-1,1:grd%jm-1)                                                       &
                                  * grd%msk(2:grd%im-1,1:grd%jm-1,k)*grd%msk(1:grd%im-2,1:grd%jm-1,k)

       grd%uvl(2:grd%im  ,1:grd%jm-1,k) = grd%uvl(2:grd%im  ,1:grd%jm-1,k) + div(1:grd%im-1,1:grd%jm-1)    &
                                          * grd%dy(2:grd%im  ,1:grd%jm-1)                                  &
                                          / grd%dx(1:grd%im-1,1:grd%jm-1) / grd%dy(1:grd%im-1,1:grd%jm-1)
       grd%uvl(1:grd%im-1,1:grd%jm-1,k) = grd%uvl(1:grd%im-1,1:grd%jm-1,k) - div(1:grd%im-1,1:grd%jm-1)    &
                                          * grd%dy(1:grd%im-1,1:grd%jm-1)                                  &
                                          / grd%dx(1:grd%im-1,1:grd%jm-1) / grd%dy(1:grd%im-1,1:grd%jm-1)
       grd%vvl(1:grd%im-1,2:grd%jm  ,k) = grd%vvl(1:grd%im-1,2:grd%jm  ,k) + div(1:grd%im-1,1:grd%jm-1)    &
                                          * grd%dx(1:grd%im-1,2:grd%jm  )                                  &
                                          / grd%dx(1:grd%im-1,1:grd%jm-1) / grd%dy(1:grd%im-1,1:grd%jm-1)
       grd%vvl(1:grd%im-1,1:grd%jm-1,k) = grd%vvl(1:grd%im-1,1:grd%jm-1,k) - div(1:grd%im-1,1:grd%jm-1)    &
                                          * grd%dx(1:grd%im-1,1:grd%jm-1)                                  &
                                          / grd%dx(1:grd%im-1,1:grd%jm-1) / grd%dy(1:grd%im-1,1:grd%jm-1)
       div(:,:) = 0.0

      enddo

     enddo

end subroutine div_dmp_ad
