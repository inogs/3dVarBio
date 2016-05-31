subroutine def_grd

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
! Define the grid                                                      !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str

  implicit none

  INTEGER(I4)    :: i, j, k
  INTEGER :: indSupWP
! ---
! Define grid 
 grd%grd_mod  = drv%grid (drv%ktr)
 grd%read_grd = drv%read_grd (drv%ktr)

!Read grid definition
      call rdgrd
     grd%dxdy(:,:) =  grd%dy(:,:) * grd%dx(:,:)

     grd%dlt = (grd%lat(1,2) - grd%lat(1,1))
     grd%dln = (grd%lon(2,1) - grd%lon(1,1))

! Define grid for horizontal covariances
     if( drv%mask(drv%ktr).eq.1)then
         grd%msr(:,:,:) = 1.0
     else if( drv%mask(drv%ktr).eq.2)then
        do k=1,grd%km
         grd%msr(:,:,k) = grd%msk(:,:,1)
        enddo

     else if( drv%mask(drv%ktr).eq.3)then
do i=1,grd%im
do j=1,grd%jm
do k=1,grd%km
grd%msr(i,j,k) = grd%msk(i,j,k)
enddo
enddo
enddo
!         grd%msr(:,:,:) = grd%msk(:,:,:)
     else
         write(drv%dia,*)'Wrong mask for horizontal covariances ',  &
                         drv%mask(drv%ktr)
         !stop
         call f_exit(21)
     endif
      grd%adxdy = sum(grd%dxdy) / (grd%im * grd%jm)

      grd%dxdy(:,:) = sqrt(grd%dxdy(:,:))
      grd%adxdy = sqrt(grd%adxdy)

      grd%nps = grd%im*grd%jm*grd%km

     do k=1,grd%km
        grd%ums(1:grd%im-1,1:grd%jm,k) =                      &
          grd%msk(1:grd%im-1,1:grd%jm,k) * grd%msk(2:grd%im,1:grd%jm,k)
        grd%vms(1:grd%im,1:grd%jm-1,k) =                      &
          grd%msk(1:grd%im,1:grd%jm-1,k) * grd%msk(1:grd%im,2:grd%jm,k)
     enddo

     do j= 1,grd%jm
      do i= 1,grd%im
       grd%f(i,j) = 0.00014584 *                &
                    sin((grd%lat(1,1)+(j-1.+0.5)*(grd%lat(2,1)-grd%lat(1,1)))*3.141592654/180.)
      enddo
     enddo


nSurfaceWaterPoints = 0
do i=1,grd%im
do j=1,grd%jm
   if (grd%msk(i,j,1).eq.1) nSurfaceWaterPoints = nSurfaceWaterPoints+1
enddo
enddo


ALLOCATE (SurfaceWaterpoints(2,nSurfaceWaterpoints))

write(*,*) 'nSurfaceWaterpoints = ', nSurfaceWaterpoints

indSupWP=0
do i=1,grd%im
do j=1,grd%jm
    if (grd%msk(i,j,1).eq.1) then
        indSupWP = indSupWP+1
        SurfaceWaterPoints(1,indSupWP) = i
        SurfaceWaterPoints(2,indSupWP) = j
    endif
enddo
enddo

end subroutine def_grd
