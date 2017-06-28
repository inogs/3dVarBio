subroutine veof_nut_ad(NutArrayAd, Var)

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

 INTEGER(i4)             :: i, j, k, l, n, k1, offset
 REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm
 REAL(r8) :: NutArrayAd(grd%im,grd%jm,grd%km)
 CHARACTER :: Var
 INTEGER   :: MyNEofs

  offset = 0
  if(Var .eq. 'N') then
      MyNEofs = ros%neof_n3n
      offset = ros%neof_chl
  else
      MyNEofs = ros%neof_o2o
      offset = ros%neof_chl + ros%neof_n3n
  endif

  do n=1,MyNEofs
    grd%ro_ad(:,:,n+offset) = 0.0 ! OMP
  enddo

!$OMP PARALLEL  &
!$OMP PRIVATE(i,j,k,k1,n) &
!$OMP PRIVATE(egm) 
!$OMP DO
  do n=1,MyNEofs

   egm(:,:) = 0.0

   ! 3D variables
   k1 = 0
   
   do k=1,grd%km ! OMP
      k1 = k1 + 1
      do j=1,grd%jm
         do i=1,grd%im
            if(Var .eq. 'N') then
                  egm(i,j) = egm(i,j) + ros%evc_n3n(grd%reg(i,j), k,n) * NutArrayAd(i,j,k)
            else
                  egm(i,j) = egm(i,j) + ros%evc_o2o(grd%reg(i,j), k,n) * NutArrayAd(i,j,k)
            endif
         enddo
      enddo
   enddo
   
   
   do j=1,grd%jm
      do i=1,grd%im
         if(Var .eq. 'N') then
            egm(i,j) = ros%eva_n3n(grd%reg(i,j),n) * egm(i,j) 
         else
            egm(i,j) = ros%eva_o2o(grd%reg(i,j),n) * egm(i,j) 
         endif         
      enddo
   enddo
   
   do j=1,grd%jm
      do i=1,grd%im
         grd%ro_ad(i,j,n+offset) = grd%ro_ad(i,j,n+offset) + egm(i,j) 
      enddo
   enddo
   
enddo
!$OMP END DO
!$OMP END PARALLEL 


end subroutine veof_nut_ad
