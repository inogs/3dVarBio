subroutine rcfl_y_ad( im, jm, km, jmax, al, bt, fld, jnx, jmx)
  
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
  !    MERCHANTABILITY or FITNESS FOR a_rcy PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a_rcy copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !--------------------------------------------------------------------------- 
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Recursive filter in y direction - adjoint
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
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
           c_rcy(i,jnx(i,j,k),tid) = fld(i,j,k)
        enddo
     enddo
     alp_rcy(:,:,tid) = al(:,:,k)
     bta_rcy(:,:,tid) = bt(:,:,k)
     
     do ktr = 1,rcf%ntr
        
        ! negative direction 
        b_rcy(:,:,tid) = 0.0
        
        do j=1,jmx(k)-1
           c_rcy(:,j+1,tid) = c_rcy(:,j+1,tid) + bta_rcy(:,j,tid)*c_rcy(:,j,tid)
           b_rcy(:,j,tid)   = (1.-bta_rcy(:,j,tid))*c_rcy(:,j,tid)
        enddo
        
        
        if( ktr.eq.1 )then
           b_rcy(:,jmx(k),tid) = b_rcy(:,jmx(k),tid) + c_rcy(:,jmx(k),tid) / (1.+bta_rcy(:,jmx(k),tid))
        else
           b_rcy(:,jmx(k),tid  ) = b_rcy(:,jmx(k),tid  ) + (1.-bta_rcy(:,jmx(k),tid)) * c_rcy(:,jmx(k),tid) / (1.-bta_rcy(:,jmx(k),tid)**2)**2
           b_rcy(:,jmx(k)-1,tid) = b_rcy(:,jmx(k)-1,tid) - (1.-bta_rcy(:,jmx(k),tid)) * &
                bta_rcy(:,jmx(k),tid)**3 * c_rcy(:,jmx(k),tid) / (1.-bta_rcy(:,jmx(k),tid)**2)**2
        endif
        
        ! positive direction 
        a_rcy(:,:,tid) = 0.0
        
        do j=jmx(k),2,-1
           b_rcy(:,j-1,tid) = b_rcy(:,j-1,tid) + alp_rcy(:,j,tid)*b_rcy(:,j,tid)
           a_rcy(:,j,tid) = a_rcy(:,j,tid) + (1.-alp_rcy(:,j,tid))*b_rcy(:,j,tid)
        enddo
        
        
        if( ktr.eq.1 )then
           a_rcy(:,1,tid) = a_rcy(:,1,tid) + (1.-alp_rcy(:,1,tid)) * b_rcy(:,1,tid)
        elseif( ktr.eq.2 )then
           a_rcy(:,1,tid) = a_rcy(:,1,tid) + b_rcy(:,1,tid) / (1.+alp_rcy(:,1,tid))
        else
           a_rcy(:,1,tid) = a_rcy(:,1,tid) + (1.-alp_rcy(:,1,tid)) * b_rcy(:,1,tid) / (1.-alp_rcy(:,1,tid)**2)**2
           a_rcy(:,2,tid) = a_rcy(:,2,tid) - (1.-alp_rcy(:,1,tid)) * alp_rcy(:,1,tid)**3 * b_rcy(:,1,tid) / (1.-alp_rcy(:,1,tid)**2)**2
        endif
        
        
        c_rcy(:,:,tid) = a_rcy(:,:,tid)
        
     enddo
     
     ! This way fills land points with some values.
     ! We prefer not investigate at the mooment and use only the water points
     do j=1,GlobalCol
        do i=1,localRow
           if(grd%global_msk(i + GlobalRowOffset,j,1).eq.1) then
              fld(i,j,k) = c_rcy(i,jnx(i,j,k),tid)
           end if
        end do
     end do

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  
end subroutine rcfl_y_ad
