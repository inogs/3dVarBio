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
  use mpi_str
  
  implicit none
  
  INTEGER(i4)     :: i, j, k, l,n, k1, my_km, MyNEofs, ierr
  REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm
  REAL(r8), ALLOCATABLE, DIMENSION(:,:)  :: eva
  REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)  :: evc
  
  my_km = 0
  if((drv%chl_assim .eq.1) .and. (drv%multiv .eq. 0)) then
    MyNEofs = ros%neof_chl
    my_km = grd%km
  else if((drv%chl_assim .eq.0) .and. (drv%multiv .eq. 1)) then
    MyNEofs = ros%neof_multi
    my_km = ros%kmchl
  endif

  if(my_km .eq. 0) then 
    if(MyId .eq. 0) then
      write(drv%dia,*) "Error! my_km for chlorophyll not setted"
      write(drv%dia,*) "chl_assim e multiv flags should be alternatively valid"
      write(drv%dia,*) ""
      write(*,*) "Error! my_km for chlorophyll not setted! Aborting"
      write(*,*) "chl_assim e multiv flags should be alternatively valid"
      write(*,*) ""
    endif
    call MPI_Barrier(Var3DCommunicator, ierr)
    call MPI_Abort(Var3DCommunicator,-1,ierr)
  endif

  ALLOCATE (eva(ros%nreg,MyNEofs)); eva = huge(eva(1,1))
  ALLOCATE (evc(ros%nreg,my_km,MyNEofs)); evc = huge(evc(1,1,1))
  
  if((drv%chl_assim .eq.1) .and. (drv%multiv .eq. 0)) then
    eva(:,:) = ros%eva_chl(:,:)
    evc(:,:,:) = ros%evc_chl(:,:,:)
  else if((drv%chl_assim .eq.0) .and. (drv%multiv .eq. 1)) then
    eva(:,:) = ros%eva_multi(:,:)
    evc(:,1:my_km,:) = ros%evc_multi(:,1:my_km,:)
  endif
  
  grd%chl(:,:,:) = 0.0
  
  !cdir noconcur
  do n=1,MyNEofs
     
     egm(:,:) = 0.0
     
     do j=1,grd%jm
        do i=1,grd%im
           egm(i,j) = eva(grd%reg(i,j),n) * grd%ro( i, j, n)
        enddo
     enddo
          
     ! 3D variables
     do k=1,my_km ! OMP
      k1 = k1 + 1
        do j=1,grd%jm
          do i=1,grd%im
            grd%chl(i,j,k) = grd%chl(i,j,k) + evc(grd%reg(i,j),k,n) * egm(i,j)
          enddo
        enddo
     enddo
  enddo

  DEALLOCATE(eva,evc)
  
end subroutine veof_chl
