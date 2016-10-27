subroutine oceanvar
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006, 2007 Srdjan Dobricic, CMCC, Bologna                   !
  !                                                                          !
  !    This file is part of OceanVar.                                        !
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
  
#ifdef _USE_MPI
  use mpi_str
#endif
  
  implicit none
  
  INTEGER(i4)   ::  ktr
  
#ifdef _USE_MPI

  INTEGER(i4)   ::  MyID
  
  call mynode
  MyID = MyRank
  
#endif
     
  ! ---
  ! Initialize diagnostics and read namelists
  call def_nml
  ! ---
  ! Outer loop - multigrid
  do ktr = 1,drv%ntr
     drv%ktr = ktr
     
     ! ---
     ! Define grid parameters
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 )then

#ifdef _USE_MPI
        call parallel_def_grd
        
        if(MyRank .eq. 0) &
#else
             
             call def_grd
#endif
        
        write(drv%dia,*) 'out of def_grd '           

     endif
     
     ! ---
     ! Get observations
     if(ktr.eq.1) call get_obs
     
#ifdef _USE_MPI
     if(MyRank .eq. 0) &
#endif
        write(drv%dia,*) 'out of get_obs'
     
     ! ---
     ! Define interpolation parameters
     call int_par
     
#ifdef _USE_MPI
     if(MyRank .eq. 0) &
#endif
          write(drv%dia,*) 'out of int_par'
             
     ! ---
     ! Define observational vector
     call obs_vec
     
#ifdef _USE_MPI
     if(MyRank .eq. 0) &
#endif
          write(drv%dia,*) 'out of obs_vec'
             
     ! ---
     ! Define constants for background covariances
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 ) then
        
#ifdef _USE_MPI
        ! call parallel_def_cov
        call def_cov
        if(MyRank .eq. 0) &
#else
        call def_cov
#endif
             write(drv%dia,*) 'out of def_cov '

     endif
     
     ! ---
     ! Initialize cost function and its gradient
     call ini_cfn

#ifdef _USE_MPI
     if(MyRank .eq. 0) &
#endif
          write(drv%dia,*) 'out of ini_cfn'

     ! ---
     ! Calculate the initial norm the gradient
     if( ktr.gt.1 ) then
        call ini_nrm
        
#ifdef _USE_MPI
        if(MyRank .eq. 0) &
#endif
             write(drv%dia,*) 'out of ini_nrm '
     endif
           
     ! ---
     ! Initialise from old iterration
     if( ktr.gt.1 .and. drv%ratio(ktr).ne.1.0 ) then
        call ini_itr
        
#ifdef _USE_MPI
        if(MyRank .eq. 0) &
#endif
             write(drv%dia,*) 'out of ini_itr '
           
     endif
     
     ! ---
     ! Minimize the cost function (inner loop)
#ifndef _USE_MPI
     call min_cfn
     write(drv%dia,*) 'out of min_cfn'
#else
     call tao_minimizer
     if(MyRank .eq. 0) then
        write(drv%dia,*) 'out of tao_minimizer'
     endif
#endif
        
     if(ktr.eq.drv%ntr)then
        ! ---
        ! Convert to innovations
        call cnv_inn
        ! ---
        ! Write outputs and diagnostics
#ifdef _USE_MPI
        ! call parallel_wrt_dia
        call wrt_dia
#else
        call wrt_dia
#endif
     endif
     
     ! ---
     ! Save old iterration
     !   if( ktr.ne.drv%ntr)then
     !    if(drv%ratio(ktr+1).ne.1.0 ) then
     call sav_itr
     
#ifdef _USE_MPI
     if(MyRank .eq. 0) &
#endif
          write(drv%dia,*) 'out of sav_itr '
        
     !    endif
     !   endif
     
     ! ---
     ! End of outer loop
     
  enddo
  !-----------------------------------------------------------------
  
  ! completely clean memory
  call clean_mem
  
  !-----------------------------------------------------------------
#ifdef _USE_MPI
  if(MyRank .eq. 0) close(drv%dia)
  call mpi_stop

#else
  close(drv%dia)
#endif
  
end subroutine oceanvar
