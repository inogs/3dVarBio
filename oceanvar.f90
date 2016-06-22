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

#ifdef key_mpp_mpi
 use myalloc_mpi
#endif

 implicit none

 INTEGER(i4)   ::  ktr

#ifdef key_mpp_mpi
 INTEGER(i4)   ::  MyID
 
 call mynode
 MyID = rank !mynode()

 print*,MyID, rank, size

 if( MyID .eq. 0) then
 print*, "Only process ", MyID, " will works..."
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
       call def_grd
          
       write(drv%dia,*) 'out of def_grd '
      endif

! ---
! Get observations
      if(ktr.eq.1) call get_obs
       write(drv%dia,*) 'out of get_obs'

! ---
! Define interpolation parameters
      call int_par
       write(drv%dia,*) 'out of int_par'

! ---
! Define observational vector
      call obs_vec
       write(drv%dia,*) 'out of obs_vec'

! ---
! Define constants for background covariances
     if( ktr.eq.1 .or. drv%ratio(ktr).ne.1.0 ) then
      call def_cov
       write(drv%dia,*) 'out of def_cov '
     endif

! ---
! Initialize cost function and its gradient
      call ini_cfn
       write(drv%dia,*) 'out of ini_cfn'

! ---
! Calculate the initial norm the gradient
     if( ktr.gt.1 ) then
      call ini_nrm
       write(drv%dia,*) 'out of ini_nrm '
     endif

! ---
! Initialise from old iterration
     if( ktr.gt.1 .and. drv%ratio(ktr).ne.1.0 ) then
      call ini_itr
       write(drv%dia,*) 'out of ini_itr '
     endif

! ---
! Minimize the cost function (inner loop)
     call min_cfn
       write(drv%dia,*) 'out of min_cfn'

   if(ktr.eq.drv%ntr)then
! ---
! Convert to innovations
     call cnv_inn
! ---
! Write outputs and diagnostics
     call wrt_dia
   endif

! ---
! Save old iterration
!   if( ktr.ne.drv%ntr)then
!    if(drv%ratio(ktr+1).ne.1.0 ) then
     call sav_itr
       write(drv%dia,*) 'out of sav_itr '
!    endif
!   endif

! ---
! End of outer loop

  enddo
!-----------------------------------------------------------------
close(drv%dia)
#ifdef key_mpp_mpi
endif
#endif

#ifdef key_mpp_mpi
CALL mpi_stop
#endif
end subroutine oceanvar
