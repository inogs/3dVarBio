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
  
#ifdef _USE_MPI
  use mpi_str
#endif
  
  implicit none
  
  INTEGER(I4)    :: i, j, k
  INTEGER :: indSupWP
  ! ---
  ! Define grid 
  grd%grd_mod  = drv%grid (drv%ktr)
  
  !Read grid definition
#ifndef _USE_MPI
  call rdgrd
#else
  ! call rdgrd
  call parallel_rdgrd
#endif
  
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
#ifdef _USE_MPI
     if(MyRank .eq. 0) then
#endif
        
     write(drv%dia,*)'Wrong mask for horizontal covariances ',  &
          drv%mask(drv%ktr)

     !stop
#ifdef _USE_MPI
  endif
  
  call MPI_Abort(MPI_COMM_WORLD, -1, i)

#else

  call f_exit(21)

#endif
  endif
  
  
  nSurfaceWaterPoints = 0
  do i=1,grd%im
     do j=1,grd%jm
        if (grd%msk(i,j,1).eq.1) nSurfaceWaterPoints = nSurfaceWaterPoints+1
     enddo
  enddo
  
  
  ALLOCATE (SurfaceWaterPoints(2,nSurfaceWaterPoints))

! #ifdef _USE_MPI
!   if(MyRank .eq. 0) &
! #endif
  write(*,*) 'nSurfaceWaterPoints = ', nSurfaceWaterPoints
  
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
