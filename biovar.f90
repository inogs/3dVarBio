subroutine biovar
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006, 2007 Srdjan Dobricic, CMCC, Bologna                   !
  !                                                                          !
  !    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                        !
  !                                                                          !
  !    3DVarBio is  free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    3DVarBio is  distributed in the hope that it will be useful,           !
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
  !                                                                          !
  !---------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! The main driver for the OceanVar                                     !
  !                                                                      !
  ! Version 0.1: S.Dobricic 2006                                         !
  !   Horizontal covariance with recursive filters, vertical with EOFs,  !
  !   assimilation of satellite observations of SLA, in situ observations!
  !   by XBT and ARGO floats                                             !
  !                                                                      !
  ! Version 0.2: S.Dobricic 2007                                         !
  !   Multigrid method. Internal boundaries for horizontal covariances.  !
  !                                                                      !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use mpi_str
  
  implicit none
  
  REAL(r8)      ::  tstart, tend
  
  call my_3dvar_node
  
  ! ---
  ! Initialize diagnostics and read namelists
  call def_nml

  ! ---
  ! Define grid parameters
  call def_grd
  if(MyId .eq. 0) write(drv%dia,*) 'out of def_grd '           
  
  ! ---
  ! Get observations
  call get_obs
  if(MyId .eq. 0) write(drv%dia,*) 'out of get_obs'
     
  ! ---
  ! Define interpolation parameters
  call int_par
  if(MyId .eq. 0) write(drv%dia,*) 'out of int_par'
             
  ! ---
  ! Define observational vector
  call obs_vec
  if(MyId .eq. 0) write(drv%dia,*) 'out of obs_vec'
             
  ! ---
  ! Define constants for background covariances
  call def_cov
  if(MyId .eq. 0) write(drv%dia,*) 'out of def_cov '
     
  ! ---
  ! Initialize cost function and its gradient
  call ini_cfn
  if(MyId .eq. 0) write(drv%dia,*) 'out of ini_cfn'

  ! ---
  ! Minimize the cost function (inner loop)
  tstart = MPI_Wtime()
  call tao_minimizer
  tend = MPI_Wtime()
  if(MyId .eq. 0) then
    write(drv%dia,*) 'out of tao_minimizer'
    write(drv%dia,*) 'minimization executed in', tend-tstart,'sec'
    write(drv%dia,*) 'time/iteration = ', (tend-tstart)/drv%MyCounter,'sec'

    write(*,*) 'minimization executed in', tend-tstart,'sec'
    write(*,*) 'time/iteration = ', (tend-tstart)/drv%MyCounter,'sec'
  endif
        
  ! ---
  ! Convert to innovations
  call cnv_inn
  ! ---
  ! Write outputs and diagnostics
  if(drv%bio_assim .eq. 0) then
    call wrt_dia
  else
    call wrt_bio_stat
  endif

  call sav_itr
  if(MyId .eq. 0) write(drv%dia,*) 'out of sav_itr '
  
  ! clean memory
  call clean_mem
  
  !-----------------------------------------------------------------
  if(MyId .eq. 0) close(drv%dia)

end subroutine biovar
