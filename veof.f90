subroutine veof
!anna
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
! Vertical transformation                           
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use grd_str
 use eof_str

 implicit none

 INTEGER(i4)     :: i, j, k, l,n, k1
 REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm


 if(drv%biol.eq.0 .or. drv%bphy.eq.1) then
     grd%eta(:,:  ) = 0.0
     grd%tem(:,:,:) = 0.0
     grd%sal(:,:,:) = 0.0
 endif
 if(drv%biol.eq.1)then
     do l=1,grd%nchl
!$OMP PARALLEL  &
!$OMP PRIVATE(k)
!$OMP DO
     do k=1,grd%km
        grd%chl(:,:,k,l) = 0.0
     enddo
!$OMP END DO
!$OMP END PARALLEL
     enddo
 endif

!cdir noconcur
  do n=1,ros%neof

     egm(:,:) = 0.0

! egm
   do l=1,n   
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
!     egm(i,j) = egm(i,j) + ros%cor( i, j, l, n) * grd%ro(i,j,l)
#else
!     egm(i,j) = egm(i,j) + ros%cor(grd%reg(i,j), l, n) * grd%ro(i,j,l)
#endif
    enddo
   enddo
   enddo


   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     egm(i,j) = ros%eva( i, j, n) *  grd%ro( i, j, n)
#else
     egm(i,j) = ros%eva(grd%reg(i,j),n) * grd%ro( i, j, n)
#endif
    enddo
   enddo


 if(drv%biol.eq.0 .or. drv%bphy.eq.1) then

! Eta
  if(drv%bmd(drv%ktr) .ne. 1)then
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     grd%eta(i,j) = grd%eta(i,j) + ros%evc( i, j,1,n) * egm(i,j)
#else
     grd%eta(i,j) = grd%eta(i,j) + ros%evc(grd%reg(i,j),1,n) * egm(i,j)
#endif
    enddo
   enddo
  endif




! 3D variables
  do k=1,grd%km
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     grd%tem(i,j,k) = grd%tem(i,j,k) + ros%evc( i, j, k+1       , n)  * egm(i,j)
     grd%sal(i,j,k) = grd%sal(i,j,k) + ros%evc( i, j, k+grd%km+1, n)  * egm(i,j)
#else
     grd%tem(i,j,k) = grd%tem(i,j,k) + ros%evc(grd%reg(i,j), k+1      ,n) * egm(i,j)
     grd%sal(i,j,k) = grd%sal(i,j,k) + ros%evc(grd%reg(i,j),k+grd%km+1,n) * egm(i,j)
#endif
    enddo
   enddo
  enddo

  write(6,*) "this branch is not parallelized!!!"
  stop
 endif

 if(drv%biol.eq.1)then

! 3D variables
   if(drv%bphy.eq.1)then
    k1 = 2*grd%km+1
   else
    k1 = 0
   endif
 do l=1,grd%nchl
!$OMP PARALLEL  &
!$OMP PRIVATE(i,j,k,k1) 
!$OMP DO
  do k=1,grd%km ! OMP
    k1 = k1 + 1
   do j=1,grd%jm
    do i=1,grd%im
#ifdef opt_huge_memory
     grd%chl(i,j,k,l) = grd%chl(i,j,k,l) + ros%evc( i, j, k1, n)  * egm(i,j)
#else
     grd%chl(i,j,k,l) = grd%chl(i,j,k,l) + ros%evc(grd%reg(i,j),k,n) * egm(i,j)
#endif
    enddo
   enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

 endif

  enddo

end subroutine veof
