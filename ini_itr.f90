subroutine ini_itr

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
  use eof_str
  use ctl_str

  implicit none

  INTEGER(i4)               :: i, j, k , ii, jj, kk
  REAL(r8)                  :: ri, rj, p, q
  REAL(r8)                  :: div_x, div_y
  REAL(r8),    ALLOCATABLE  :: pq1(:,:), pq2(:,:), pq3(:,:), pq4(:,:)
  INTEGER(i4), ALLOCATABLE  :: i1(:,:), j1(:,:)



   ALLOCATE( pq1(grd%im,grd%jm) );  pq1 = huge(pq1(1,1))
   ALLOCATE( pq2(grd%im,grd%jm) );  pq2 = huge(pq2(1,1))
   ALLOCATE( pq3(grd%im,grd%jm) );  pq3 = huge(pq3(1,1))
   ALLOCATE( pq4(grd%im,grd%jm) );  pq4 = huge(pq4(1,1))
   ALLOCATE(  i1(grd%im,grd%jm) );  i1  = huge(i1 (1,1))
   ALLOCATE(  j1(grd%im,grd%jm) );  j1  = huge(j1 (1,1))

! ---
! Interpolate between grids

     do jj=1,grd%jm
     do ii=1,grd%im
       ri=max(1.,min(real(drv%im-1),real(ii-1)/real(drv%ratio(drv%ktr)) + 1.))
        i=int(ri)
       p=ri-i
       rj=max(1.,min(real(drv%jm-1),real(jj-1)/real(drv%ratio(drv%ktr)) + 1.))
        j=int(rj)
       q=rj-j

           i1(ii,jj) = i
           j1(ii,jj) = j

         div_y =  (1.-q) * max(drv%msk(i,j  ),drv%msk(i+1,j  ))      &
                 +    q  * max(drv%msk(i,j+1),drv%msk(i+1,j+1))
         div_x =  (1.-p) * drv%msk(i  ,j) + p * drv%msk(i+1,j)
          pq1(ii,jj) = drv%msk(i,j)                                  &
                      * max(drv%msk(i,j),drv%msk(i+1,j))             &
                       * (1.-p) * (1.-q)                             &
                      /( div_x * div_y + 1.e-16 )
          pq2(ii,jj) = drv%msk(i+1,j)                                &
                      * max(drv%msk(i,j),drv%msk(i+1,j))             &
                      *     p  * (1.-q)                              &
                      /( div_x * div_y + 1.e-16 )
         div_x =  (1.-p) * drv%msk(i  ,j+1) + p * drv%msk(i+1,j+1)
          pq3(ii,jj) = drv%msk(i,j+1)                                &
                      * max(drv%msk(i,j+1),drv%msk(i+1,j+1))         &
                      * (1.-p) *     q                               &
                      /( div_x * div_y + 1.e-16 )
          pq4(ii,jj) = drv%msk(i+1,j+1)                              &
                      * max(drv%msk(i,j+1),drv%msk(i+1,j+1))         &
                      *     p  *     q                               &
                      /( div_x * div_y + 1.e-16 )

     enddo
     enddo


    do k = 1,ros%neof
     do jj=1,grd%jm
     do ii=1,grd%im
        i=i1(ii,jj)
        j=j1(ii,jj)
       grd%ro(ii,jj,k) = pq1(ii,jj) * drv%ro(i,j  ,k) + pq2(ii,jj) * drv%ro(i+1,j  ,k)  &
                       + pq3(ii,jj) * drv%ro(i,j+1,k) + pq4(ii,jj) * drv%ro(i+1,j+1,k) 
     enddo
     enddo
    enddo

! ---
! Reconstruct the control vector
       kk = 0
   do k=1,ros%neof
    do j=1,grd%jm
     do i=1,grd%im
       kk = kk+1
       ctl%x_c(kk) = grd%ro(i,j,k)/drv%ratio(drv%ktr)
     enddo
    enddo
   enddo


   DEALLOCATE( drv%ro, drv%ro_ad, drv%msk)
   DEALLOCATE( pq1 )
   DEALLOCATE( pq2 )
   DEALLOCATE( pq3 )
   DEALLOCATE( pq4 )
   DEALLOCATE(  i1 )
   DEALLOCATE(  j1 )

! Calculate the cost function and its gradient

          call costf



end subroutine ini_itr
