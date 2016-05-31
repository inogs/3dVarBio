subroutine ini_nrm


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
! Minimise the cost function                                           !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str

 implicit none

 INTEGER(i4)           :: k
 REAL(r8)              :: maxpg
 REAL(r8), allocatable, dimension (:) :: x_s, g_s

! Save in temporary arrays
 ALLOCATE ( x_s(ctl%n)) ; x_s(:) = ctl%x_c(:)
 ALLOCATE ( g_s(ctl%n)) ; g_s(:) = ctl%g_c(:)

! Initialize the control vector to zero
          ctl%x_c(:) =  0.0

! Calculate the cost function and its gradient
          call costf

! Calculate the norm of the gradient
            maxpg = 0.0
           do k=1,ctl%n
            maxpg = max(maxpg,abs(ctl%g_c(k)))
           enddo

            ctl%pgtol = maxpg * ctl%pgper

         print*,' maxpg  is: ',maxpg, ctl%pgtol

          ctl%x_c(:) = x_s(:)
          ctl%g_c(:) = g_s(:)

 DEALLOCATE ( x_s, g_s )

end subroutine ini_nrm
