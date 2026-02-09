subroutine veof_dnc
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
  use mpi_str
  
  implicit none
  
  INTEGER(i4)     :: i, j, k, l,n, my_km, MyNEofs, ierr
  REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm
  REAL(r8), ALLOCATABLE, DIMENSION(:,:)  :: eva
  REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)  :: evc
  
  my_km = grd%km
  MyNEofs = ros%neof_dnc
  offset = ros%neof_chl + ros%neof_n3n


  ALLOCATE (eva(ros%nreg,MyNEofs)); eva = huge(eva(1,1))
  ALLOCATE (evc(ros%nreg,my_km,MyNEofs)); evc = huge(evc(1,1,1))
  
  eva(:,:) = ros%eva_dnc(:,:)
  evc(:,1:my_km,:) = ros%evc_dnc(:,my_km+1:my_km*2,:)
  
  grd%dnc(:,:,:) = 0.0
  
  !cdir noconcur
  do n=1,MyNEofs
     
     egm(:,:) = 0.0
     
     do j=1,grd%jm
        do i=1,grd%im
           egm(i,j) = eva(grd%reg(i,j),n) * grd%ro( i, j, n+offset)
        enddo
     enddo
          
     ! 3D variables
     do k=1,my_km ! OMP
        do j=1,grd%jm
          do i=1,grd%im
            grd%dnc(i,j,k) = grd%dnc(i,j,k) + evc(grd%reg(i,j),k,n) * egm(i,j)
          enddo
        enddo
     enddo
  enddo

  DEALLOCATE(eva,evc)
  
end subroutine veof_dnc
