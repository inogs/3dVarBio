subroutine veof_nut(NutArray, Var)
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
  
  INTEGER(i4)     :: i, j, k, l,n, my_km, ierr
  REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm
  REAL(r8) :: NutArray(grd%im,grd%jm,grd%km)
  REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: eva
  REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: evc
  INTEGER(I4) :: MyNEofs, offset  
  CHARACTER   :: Var
  
  NutArray(:,:,:) = 0.0

  my_km = 0
  offset = 0
  if((drv%nut .eq.1) .and. (drv%multiv .eq. 0)) then
    my_km = grd%km
    if(Var .eq. 'N') then
      MyNEofs = ros%neof_n3n
      offset = ros%neof_chl
    else
      MyNEofs = ros%neof_o2o
      offset = ros%neof_chl + ros%neof_n3n
    endif
  else if((drv%nut .eq.0) .and. (drv%multiv .eq. 1)) then
    if(Var .eq. 'N') then
      my_km = grd%km
      MyNEofs = ros%neof_multi
    else
      if(MyId .eq. 0) then
        write(drv%dia,*) "Error! Only nitrate multivariate assimilation implemented"
        write(drv%dia,*) ""
        write(*,*) "Error! Only nitrate multivariate assimilation implemented! Aborting"
        write(*,*) ""
      endif
      call MPI_Barrier(Var3DCommunicator, ierr)
      call MPI_Abort(Var3DCommunicator,-1,ierr)  
    endif
  endif


  if(my_km .eq. 0) then
    if(MyId .eq. 0) then
      write(drv%dia,*) "Error! my_km for nutrient not setted"
      write(drv%dia,*) "drv%nut e multiv flags should be alternatively valid"
      write(drv%dia,*) ""
      write(*,*) "Error! my_km for nutrient not setted! Aborting"
      write(*,*) "drv%nut e multiv flags should be alternatively valid"
      write(*,*) ""
    endif
    call MPI_Barrier(Var3DCommunicator, ierr)
    call MPI_Abort(Var3DCommunicator,-1,ierr)
  endif   

  ALLOCATE (eva(ros%nreg,MyNEofs)); eva = huge(eva(1,1))
  ALLOCATE (evc(ros%nreg,my_km,MyNEofs)); evc = huge(evc(1,1,1))

  if((drv%nut .eq.1) .and. (drv%multiv .eq. 0)) then
    if(Var .eq. 'N') then
      eva = ros%eva_n3n
      evc = ros%evc_n3n
    else
      eva = ros%eva_o2o
      evc = ros%evc_o2o
    endif
  else if((drv%nut .eq.0) .and. (drv%multiv .eq. 1)) then
    if(Var .eq. 'N') then
      eva = ros%eva_multi
      evc(:,1:my_km,:) = ros%evc_multi(:,ros%kmchl+1:ros%kmchl+grd%km,:)
    endif
  endif
  
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
              NutArray(i,j,k) = NutArray(i,j,k) + evc(grd%reg(i,j),k,n) * egm(i,j)
           enddo
        enddo
     enddo
  enddo

DEALLOCATE(eva,evc)
end subroutine veof_nut
