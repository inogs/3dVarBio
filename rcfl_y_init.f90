subroutine rcfl_y_init( im, jm, km, jmax, al, bt, fld, jnx, jmx)

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
! Recursive filter in y direction
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use cns_str

 implicit none 

 INTEGER(i4)    :: im, jm, km, jmax

 REAL(r8)       :: fld(im,jm,km)
 REAL(r8)       :: al(im,jmax,km), bt(im,jmax,km)
 INTEGER(i4)    :: jnx(im,jm,km), jmx(km)
 REAL(r8), allocatable       :: a(:,:), b(:,:), c(:,:)
 REAL(r8), allocatable       :: alp(:,:), bta(:,:)

 INTEGER(i4)    :: i,j,k, ktr


   ALLOCATE ( a(im,jmax), b(im,jmax), c(im,jmax) )
   ALLOCATE ( alp(im,jmax), bta(im,jmax) ) 
   alp = huge(alp(1,1)); bta = huge(bta(1,1))
   do k=1,km

        a(:,:) = 0.0
        b(:,:) = 0.0
        c(:,:) = 0.0
  
        do j=1,jm
         do i=1,im
            a(i,jnx(i,j,k)) = fld(i,j,k)
         enddo
        enddo
         alp(:,:) = al(:,:,k)
         bta(:,:) = bt(:,:,k)


       do ktr = 1,rcf%ntr  !numero di cicli del filtro

! positive direction
         if( ktr.eq.1 )then
           b(:,1) = (1.-alp(:,1)) * a(:,1)
         elseif( ktr.eq.2 )then
            b(:,1) = a(:,1) / (1.+alp(:,1))
         else
            b(:,1) = (1.-alp(:,1)) * (a(:,1)-alp(:,1)**3 * a(:,2)) / (1.-alp(:,1)**2)**2
         endif

        do j=2,jmx(k)
             b(:,j) = alp(:,j)*b(:,j-1) + (1.-alp(:,j))*a(:,j)
        enddo

! negative direction
         if( ktr.eq.1 )then
           c(:,jmx(k)) = b(:,jmx(k)) / (1.+bta(:,jmx(k)))
         else
           c(:,jmx(k)) = (1.-bta(:,jmx(k))) * (b(:,jmx(k))-bta(:,jmx(k))**3 * b(:,jmx(k)-1)) / (1.-bta(:,jmx(k))**2)**2
         endif

         do j=jmx(k)-1,1,-1
          c(:,j) = bta(:,j)*c(:,j+1) + (1.-bta(:,j))*b(:,j)
         enddo

         a(:,:) = c(:,:)

       enddo

        do j=1,jm
         do i=1,im
          fld(i,j,k) = a(i,jnx(i,j,k)) 
         enddo
        enddo

   enddo

   DEALLOCATE ( a, b, c )
   DEALLOCATE ( alp, bta )

end subroutine rcfl_y_init
