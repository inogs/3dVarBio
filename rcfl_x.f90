subroutine rcfl_x( im, jm, km, imax, al, bt, fld, inx, imx)
  
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
  !    MERCHANTABILITY or FITNESS FOR a_rcx PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a_rcx copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !---------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Recursive filter in x direction                                      !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use cns_str
  use rcfl
  use grd_str
#ifdef _USE_MPI
  use mpi_str
#endif
  implicit none
  
  INTEGER(i4)    :: im, jm, km, imax
  
  REAL(r8)       :: fld(im,jm,km)
  REAL(r8)       :: al(jm,imax,km), bt(jm,imax,km)
  INTEGER(i4)    :: inx(im,jm,km), imx(km)
  
  
  INTEGER(i4)    :: i,j,k, ktr
  INTEGER(i4)    :: indSupWP
  INTEGER nthreads, tid
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  
  tid=1
  !$OMP PARALLEL &
  !$OMP PRIVATE(k,j,i,ktr,indSupWP,tid)
  !$ tid      = OMP_GET_THREAD_NUM()+1
  !$OMP DO
  do k=1,km
     
     a_rcx(:,:,tid) = 0.0
     b_rcx(:,:,tid) = 0.0
     c_rcx(:,:,tid) = 0.0
     
     do j=1,jm
        do i=1,im
           a_rcx(j,inx(i,j,k),tid) = fld(i,j,k)
        enddo
     enddo
     alp_rcx(:,:,tid) = al(:,:,k)
     bta_rcx(:,:,tid) = bt(:,:,k)
     
     do ktr = 1,rcf%ntr
        
        ! positive direction
        if( ktr.eq.1 )then
           b_rcx(:,1,tid) = (1.-alp_rcx(:,1,tid)) * a_rcx(:,1,tid)
        elseif( ktr.eq.2 )then
           b_rcx(:,1,tid) = a_rcx(:,1,tid) / (1.+alp_rcx(:,1,tid))
        else
           b_rcx(:,1,tid) = (1.-alp_rcx(:,1,tid)) * (a_rcx(:,1,tid)-alp_rcx(:,1,tid)**3 * a_rcx(:,2,tid)) / (1.-alp_rcx(:,1,tid)**2)**2
        endif
        
        do j=2,imx(k)
           b_rcx(:,j,tid) = alp_rcx(:,j,tid)*b_rcx(:,j-1,tid) + (1.-alp_rcx(:,j,tid))*a_rcx(:,j,tid)
        enddo
        
        ! negative direction
        if( ktr.eq.1 )then
           c_rcx(:,imx(k),tid) = b_rcx(:,imx(k),tid) / (1.+bta_rcx(:,imx(k),tid))
        else
           c_rcx(:,imx(k),tid) = (1.-bta_rcx(:,imx(k),tid)) * &
                (b_rcx(:,imx(k),tid)-bta_rcx(:,imx(k),tid)**3 * b_rcx(:,imx(k)-1,tid)) / (1.-bta_rcx(:,imx(k),tid)**2)**2
        endif
        
        do j=imx(k)-1,1,-1
           c_rcx(:,j,tid) = bta_rcx(:,j,tid)*c_rcx(:,j+1,tid) + (1.-bta_rcx(:,j,tid))*b_rcx(:,j,tid)
        enddo
        
        a_rcx(:,:,tid) = c_rcx(:,:,tid)
        
     enddo
     
     !        do j=1,jm
     !        do i=1,im
     !         fld(i,j,k) = a_rcx(j,inx(i,j,k))
     !        enddo
     !        enddo
     ! This way fills land points with some values.
     ! We prefer not investigate at the mooment and use only the water points
#ifdef _USE_MPI
!     do j=1,localCol
!        do i=1,GlobalRow
!           ! if(grd%global_msk(i,j + MyRank*localCol,1).eq.1) then
!           if(grd%global_msk(i,j + GlobalColOffset,1).eq.1) then
!              fld(i,j,k) = a_rcx(j,inx(i,j,k),tid)
!           end if
!        end do
!     end do
     do indSupWP=1,nSurfaceWaterPoints
        i = SurfaceWaterPoints(1,indSupWP)
        j = SurfaceWaterPoints(2,indSupWP)
        fld(i,j,k) = a_rcx(j,inx(i,j,k),tid)
     enddo
#else
     do indSupWP=1,nSurfaceWaterPoints
        i = SurfaceWaterPoints(1,indSupWP)
        j = SurfaceWaterPoints(2,indSupWP)
        fld(i,j,k) = a_rcx(j,inx(i,j,k),tid)
     enddo
#endif
     
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  
  

end subroutine rcfl_x
