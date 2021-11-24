subroutine veof_multiv_ad

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
! Vertical transformation (adjoint)                                    !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use grd_str
 use eof_str
 use mpi_str

 implicit none

 INTEGER(i4)             :: i, j, k, l, n, my_km !, k1
 REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm

my_km = ros%kmchl


do n=1,ros%neof_multi
  grd%ro_ad(:,:,n) = 0.0 ! OMP
enddo

!$OMP PARALLEL  &
!$OMP PRIVATE(i,j,k,k1,n) &
!$OMP PRIVATE(egm) 
!$OMP DO
do n=1,ros%neof_multi

  egm(:,:) = 0.0

   ! 3D variables
   ! k1 = 0
   
  do k=1,grd%km ! OMP
      ! k1 = k1 + 1
    do j=1,grd%jm
      do i=1,grd%im
           if(k.le.my_km) then
             egm(i,j) = egm(i,j) + ros%evc_multi(grd%reg(i,j), k,n) * grd%chl_ad(i,j,k)
  
           endif 
           egm(i,j) = egm(i,j) + ros%evc_multi(grd%reg(i,j),k+my_km,n) * grd%n3n_ad(i,j,k)
      enddo
    enddo
  enddo


   
  do j=1,grd%jm
    do i=1,grd%im
         egm(i,j) = ros%eva_multi(grd%reg(i,j),n) * egm(i,j) 
    enddo
  enddo
   
   !cdir serial
   ! 3D variables
   !  do l=n,ros%neof
  do j=1,grd%jm
    do i=1,grd%im
        grd%ro_ad(i,j,n) = grd%ro_ad(i,j,n) + egm(i,j)
    enddo
  enddo
   !  enddo
   !cdir end serial
   
enddo
!$OMP END DO
!$OMP END PARALLEL


end subroutine veof_multiv_ad
