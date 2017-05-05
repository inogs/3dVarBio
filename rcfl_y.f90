subroutine rcfl_y( im, jm, km, jmax, al, bt, fld, jnx, jmx)
  
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
  !    MERCHANTABILITY or FITNESS FOR a_rcy PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a_rcy copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !--------------------------------------------------------------------------- 
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Recursive filter in y direction
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use cns_str
  use rcfl
  use grd_str
  use mpi_str

  implicit none 
  
  INTEGER(i4)    :: im, jm, km, jmax
  
  REAL(r8)       :: fld(im,jm,km)
  REAL(r8)       :: al(im,jmax,km), bt(im,jmax,km)
  INTEGER(i4)    :: jnx(im,jm,km), jmx(km)
  
  INTEGER(i4)    :: i,j,k, ktr
  INTEGER(i4)    :: indSupWP
  INTEGER nthreads, tid
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  
  tid=1
  !$OMP PARALLEL  &
  !$OMP PRIVATE(k,j,i,ktr,indSupWP,tid)
  !$  tid      = OMP_GET_THREAD_NUM()+1
  
  !$OMP DO
  do k=1,km
     
     a_rcy(:,:,tid) = 0.0
     b_rcy(:,:,tid) = 0.0
     c_rcy(:,:,tid) = 0.0
     
     do j=1,jm
        do i=1,im
           a_rcy(i,jnx(i,j,k),tid) = fld(i,j,k)
        enddo
     enddo
     alp_rcy(:,:,tid) = al(:,:,k)
     bta_rcy(:,:,tid) = bt(:,:,k)
     
     
     do ktr = 1,rcf%ntr
        
        ! positive direction
        if( ktr.eq.1 )then
           b_rcy(:,1,tid) = (1.-alp_rcy(:,1,tid)) * a_rcy(:,1,tid)
        elseif( ktr.eq.2 )then
           b_rcy(:,1,tid) = a_rcy(:,1,tid) / (1.+alp_rcy(:,1,tid))
        else
           b_rcy(:,1,tid) = (1.-alp_rcy(:,1,tid)) * (a_rcy(:,1,tid)-alp_rcy(:,1,tid)**3 * a_rcy(:,2,tid)) / (1.-alp_rcy(:,1,tid)**2)**2
        endif
        
        do j=2,jmx(k)
           b_rcy(:,j,tid) = alp_rcy(:,j,tid)*b_rcy(:,j-1,tid) + (1.-alp_rcy(:,j,tid))*a_rcy(:,j,tid)
        enddo
        
        ! negative direction
        if( ktr.eq.1 )then
           c_rcy(:,jmx(k),tid) = b_rcy(:,jmx(k),tid) / (1.+bta_rcy(:,jmx(k),tid))
        else
           c_rcy(:,jmx(k),tid) = (1.-bta_rcy(:,jmx(k),tid)) * &
                (b_rcy(:,jmx(k),tid)-bta_rcy(:,jmx(k),tid)**3 * b_rcy(:,jmx(k)-1,tid)) / (1.-bta_rcy(:,jmx(k),tid)**2)**2
        endif
        
        do j=jmx(k)-1,1,-1
           c_rcy(:,j,tid) = bta_rcy(:,j,tid)*c_rcy(:,j+1,tid) + (1.-bta_rcy(:,j,tid))*b_rcy(:,j,tid)
        enddo

        a_rcy(:,:,tid) = c_rcy(:,:,tid)
        
     enddo
     
     !        do j=1,jm
     !         do i=1,im
     !          fld(i,j,k) = a_rcy(i,jnx(i,j,k))
     !         enddo
     !        enddo
     
     ! This way fills land points with some values.
     ! We prefer not investigate at the mooment and use only the water points

     do j=1,GlobalCol
        do i=1,localRow
           if(grd%global_msk(i + GlobalRowOffset,j,1).eq.1) then
              fld(i,j,k) = a_rcy(i,jnx(i,j,k),tid)
           end if
        end do
     end do
     
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  
end subroutine rcfl_y
