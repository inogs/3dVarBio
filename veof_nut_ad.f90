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

 INTEGER(i4)             :: i, j, k, l, n, offset, my_km
 REAL(r8), DIMENSION ( grd%im, grd%jm)  :: egm
 REAL(r8) :: NutArrayAd(grd%im,grd%jm,grd%km)
 REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: eva
 REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: evc
 CHARACTER :: Var
 INTEGER   :: MyNEofs

  my_km = 0
  ! Altrove usato grd%km come limite per assimilazione nit qui ro%kmnit
  ! Da correggere o fare un check
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
      my_km = ros%kmnit
      MyNEofs = ros%neof_multi
    endif
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
      evc(:,1:my_km,:) = ros%evc_multi(:,ros%kmchl+1:ros%kmchl+ros%kmnit,:)
    endif
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
   
   do k=1,my_km ! OMP
      do j=1,grd%jm
         do i=1,grd%im
            egm(i,j) = egm(i,j) + evc(grd%reg(i,j), k,n) * NutArrayAd(i,j,k)
         enddo
      enddo
   enddo
   
   
   do j=1,grd%jm
      do i=1,grd%im
        egm(i,j) = eva(grd%reg(i,j),n) * egm(i,j) 
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

DEALLOCATE(eva,evc)

end subroutine veof_nut_ad
