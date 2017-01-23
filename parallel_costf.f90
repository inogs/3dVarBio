subroutine parallel_costf
  
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
  ! Calclate the cost function and its gradient                          !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  use mpi_str
  
  implicit none
  integer :: ierr
  ! -------------------------------------------------------
  ! calculate backgorund cost term
  ! -------------------------------------------------------

  ctl%f_b = 0.5 * dot_product( ctl%x_c, ctl%x_c)
  call MPI_Allreduce(MPI_IN_PLACE, ctl%f_b, 1, MPI_REAL8, MPI_SUM, MyCommWorld, ierr)
  !    write(*,*) 'COSTF f_b = ', ctl%f_b
  
  ! -------------------------------------------------------
  ! calculate observational cost term
  ! -------------------------------------------------------

  ! --------
  ! Convert the control vector to v
  call cnv_ctv
  
  ! --------
  ! Control to physical space 
  call parallel_ver_hor
  
  ! --------
  ! Apply observational operators
  call obsop
  
  ! --------
  ! Calculate residuals
  call resid
  
  ! --------
  ! calculate cost
  ctl%f_o = 0.5 * dot_product( obs%amo, obs%amo)
  call MPI_Allreduce(MPI_IN_PLACE, ctl%f_o, 1, MPI_REAL8, MPI_SUM, MyCommWorld, ierr)
  
  ! -------------------------------------------------------
  ! Cost function
  ! -------------------------------------------------------
   
  ctl%f_c = ctl%f_b + ctl%f_o
  
  if(MyRank .eq. 0) &
       print*,' Cost function ',ctl%f_c, '(iter.',drv%MyCounter,')'
  
  ! -------------------------------------------------------
  ! calculate the cost function gradient
  ! -------------------------------------------------------
  
  ! --------
  ! Reset the increments
  call res_inc
  
  ! --------
  ! Observational operators
  call obsop_ad
  
  ! --------
  ! Control to physical space 
  call parallel_ver_hor_ad
  
  !   write(*,*) 'COSTF sum(ro_ad) = ' , sum(grd%ro_ad)
  ! --------
  ! Convert the control vector 
  call cnv_ctv_ad
  
  ! -------------------------------------------------------
  ! Cost function gradient
  ! -------------------------------------------------------
  
  !   write(*,*) 'COSTF sum(g_c) = ' , sum( ctl%g_c)
  ctl%g_c(:) = ctl%x_c(:) + ctl%g_c(:) ! OMP
  
end subroutine parallel_costf
