subroutine veof_ad

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

 implicit none

 INTEGER(i4)             :: i, j, k, l, n, k1
 REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm

  grd%ro_ad(:,:,:) = 0.0 ! OMP

!$OMP PARALLEL  &
!$OMP PRIVATE(i,j,k,k1,n) &
!$OMP PRIVATE(egm) 
!$OMP DO
  do n=1,ros%neof

   egm(:,:) = 0.0

 if(drv%biol.eq.0 .or. drv%bphy.eq.1) then

! Eta
  if(drv%bmd(drv%ktr) .ne. 1)then
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = egm(i,j) + ros%evc( i, j, 1,n) * grd%eta_ad(i,j)
#else
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), 1,n) * grd%eta_ad(i,j)
#endif
    enddo
   enddo
  endif



! 3D variables
  do k=1,grd%km
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = egm(i,j) + ros%evc( i, j, k+1       ,n) * grd%tem_ad(i,j,k)
     egm(i,j) = egm(i,j) + ros%evc( i, j, k+grd%km+1,n) * grd%sal_ad(i,j,k)
#else
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), k+1       ,n) * grd%tem_ad(i,j,k)
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), k+grd%km+1,n) * grd%sal_ad(i,j,k)
#endif
    enddo
   enddo
  enddo

 endif


 if(drv%biol.eq.1)then

! 3D variables
   if(drv%bphy.eq.1)then
    k1 = 2*grd%km+1
   else
    k1 = 0
   endif
 do l=1,grd%nchl
  do k=1,grd%km ! OMP
    k1 = k1 + 1
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = egm(i,j) + ros%evc( i, j, k1,n) * grd%chl_ad(i,j,k,l)
#else
     egm(i,j) = egm(i,j) + ros%evc(grd%reg(i,j), k,n) * grd%chl_ad(i,j,k,l)
#endif
    enddo
   enddo
  enddo
 enddo

 endif


   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = ros%eva( i, j, n) * egm(i,j) 
#else
     egm(i,j) = ros%eva(grd%reg(i,j),n) * egm(i,j) 
#endif
    enddo
   enddo

!cdir serial
! 3D variables
!  do l=n,ros%neof
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     grd%ro_ad(i,j,n) = grd%ro_ad(i,j,n) + egm(i,j) ! * ros%cor( i, j, n, l) 
#else
     grd%ro_ad(i,j,n) = grd%ro_ad(i,j,n) + egm(i,j) ! * ros%cor( grd%reg(i,j), n, l) 
#endif
    enddo
   enddo
!  enddo
!cdir end serial

  enddo
!$OMP END DO
!$OMP END PARALLEL 


end subroutine veof_ad
