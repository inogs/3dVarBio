subroutine veof_chl
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
  
  
  grd%chl(:,:,:) = 0.0
  
  !cdir noconcur
  do n=1,ros%neof_chl
     
     egm(:,:) = 0.0
     
     do j=1,grd%jm
        do i=1,grd%im
           egm(i,j) = ros%eva_chl(grd%reg(i,j),n) * grd%ro( i, j, n)
        enddo
     enddo
          
     ! 3D variables
     do k=1,grd%km ! OMP
      k1 = k1 + 1
        do j=1,grd%jm
          do i=1,grd%im
            grd%chl(i,j,k) = grd%chl(i,j,k) + ros%evc_chl(grd%reg(i,j),k,n) * egm(i,j)
          enddo
        enddo
     enddo
  enddo
  
end subroutine veof_chl
