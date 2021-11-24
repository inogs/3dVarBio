subroutine def_nml

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
! Define analysis parameters from namelists                            !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

  use set_knd
  use drv_str
  use grd_str
  use obs_str
  use eof_str
  use cns_str
  use ctl_str
  use mpi_str
  use bio_str

  implicit none

  LOGICAL       :: read_eof
  INTEGER(i4)   :: neof_chl, neof_n3n, neof_o2o, nreg, rcf_ntr
  INTEGER(i4)   :: neof_multi, kmchl, kmnit
  INTEGER(i4)   :: verbose
  REAL(r8)      :: rcf_L, ctl_tol, ctl_per, rcf_efc
  
  NAMELIST /ctllst/ ctl_tol, ctl_per, verbose
  NAMELIST /covlst/ neof_chl, neof_n3n, neof_o2o, neof_multi, kmchl, kmnit, nreg, read_eof, rcf_ntr, rcf_L, rcf_efc


! -------------------------------------------------------------------
! Open a formatted file for the diagnostics
! ---

  drv%dia = 12

  if(MyId .eq. 0) then
    open ( drv%dia, file='BioVar.diagnostics', form='formatted' )
  endif

!---------------------------------------------------------------------
! Open the namelist
! ---

  open(11,file='var_3d_nml',form='formatted')

  ! ---
  read(11,ctllst)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '  '
    write(drv%dia,*) '                      NAMELISTS: '
    write(drv%dia,*) '  '
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' MINIMIZER NAMELIST INPUT: '
    write(drv%dia,*) ' Minimum gradient of J:           ctl_tol  = ', ctl_tol
    write(drv%dia,*) ' Percentage of initial gradient:  ctl_per  = ', ctl_per
    write(drv%dia,*) ' Add verbose on standard output   verbose  = ', verbose

  endif

  ctl%pgtol = ctl_tol
  ctl%pgper = ctl_per
  drv%Verbose = verbose

! ---
  read(11,covlst)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' COVARIANCE NAMELIST INPUT: '
    write(drv%dia,*) ' Number of EOFs for chl:          neof_chl = ', neof_chl
    write(drv%dia,*) ' Number of EOFs for N3n:          neof_n3n = ', neof_n3n
    write(drv%dia,*) ' Number of EOFs for O2o:          neof_o2o = ', neof_o2o
    write(drv%dia,*) ' Number of multivariate EOFs:   neof_multi = ', neof_multi
    write(drv%dia,*) ' Chl Nlevels in multi EOFs:          kmchl = ', kmchl
    write(drv%dia,*) ' Nit Nlevels in multi EOFs:          kmnit = ', kmnit
    write(drv%dia,*) ' Number of regions:               nreg     = ', nreg
    write(drv%dia,*) ' Read EOFs from a file:           read_eof = ', read_eof
    write(drv%dia,*) ' Half number of iterations:       rcf_ntr  = ', rcf_ntr
    write(drv%dia,*) ' Horizontal correlation radius:   rcf_L    = ', rcf_L
    write(drv%dia,*) ' Extension factor for coastlines: rcf_efc  = ', rcf_efc

  endif

  close(11)

  ros%neof_chl = neof_chl
  ros%neof_n3n = neof_n3n
  ros%neof_o2o = neof_o2o
  ros%neof_multi = neof_multi
  ros%kmchl    = kmchl
  ros%kmnit    = kmnit
  ros%nreg     = nreg
  ros%read_eof = read_eof
  rcf%ntr      = rcf_ntr
  rcf%L        = rcf_L
  rcf%efc      = rcf_efc


end subroutine def_nml
