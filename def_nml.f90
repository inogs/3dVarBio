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

  LOGICAL       :: read_eof, ApplyConditions
  INTEGER(i4)   :: neof, nreg, rcf_ntr
  INTEGER(i4)   :: ctl_m, bio_assim
  INTEGER(i4)   :: biol, bphy, nphyto, uniformL, anisL, verbose
  REAL(r8)      :: rcf_L, ctl_tol, ctl_per, rcf_efc, chl_dep
  INTEGER(i4)   :: argo, sat_obs, ncmp
  
  NAMELIST /ctllst/ ctl_tol, ctl_per
  NAMELIST /covlst/ neof, nreg, read_eof, rcf_ntr, rcf_L, rcf_efc
  NAMELIST /biolst/ bio_assim, nphyto, chl_dep, ncmp, ApplyConditions
  NAMELIST /params/ sat_obs, argo, uniformL, anisL, verbose


! -------------------------------------------------------------------
! Open a formatted file for the diagnostics
! ---

  drv%dia = 12

  if(MyId .eq. 0) then
    open ( drv%dia, file='OceanVar.diagnostics', form='formatted' )
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
    write(drv%dia,*) ' Number of saved vectors:         ctl_m    = ', ctl_m
    write(drv%dia,*) ' Minimum gradient of J:           ctl_tol  = ', ctl_tol
    write(drv%dia,*) ' Percentage of initial gradient:  ctl_per  = ', ctl_per

  endif

  ! ctl%m     = ctl_m
  ctl%pgtol = ctl_tol
  ctl%pgper = ctl_per

! ---
  read(11,covlst)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' COVARIANCE NAMELIST INPUT: '
    write(drv%dia,*) ' Number of EOFs:                  neof     = ', neof
    write(drv%dia,*) ' Number of regions:               nreg     = ', nreg
    write(drv%dia,*) ' Read EOFs from a file:           read_eof = ', read_eof
    write(drv%dia,*) ' Half number of iterations:       rcf_ntr  = ', rcf_ntr
    write(drv%dia,*) ' Horizontal correlation radius:   rcf_L    = ', rcf_L
    write(drv%dia,*) ' Extension factor for coastlines: rcf_efc  = ', rcf_efc

  endif

  ros%neof     = neof
  ros%nreg     = nreg
  ros%read_eof = read_eof
  rcf%ntr      = rcf_ntr
  rcf%L        = rcf_L
  rcf%efc      = rcf_efc

! ---
  read(11,biolst)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' BIOLOGY NAMELIST INPUT: '
    write(drv%dia,*) ' Biological repartition of the chlorophyll = ', bio_assim
    write(drv%dia,*) ' Number of phytoplankton species  nphyt    = ', nphyto
    write(drv%dia,*) ' Minimum depth for chlorophyll    chl_dep  = ', chl_dep
    write(drv%dia,*) ' Number of phytoplankton components  ncmp  = ', ncmp
    write(drv%dia,*) ' Apply conditions flag     ApplyConditions = ', ApplyConditions

  endif

  drv%bio_assim = bio_assim
  bio%nphy = nphyto
  sat%dep  = chl_dep
  bio%ncmp = ncmp
  bio%ApplyConditions = ApplyConditions

  read(11,params)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' PARAMETERS NAMELIST INPUT: '
    write(drv%dia,*) ' Read Satellite observations      sat_obs  = ', sat_obs
    write(drv%dia,*) ' Read ARGO float observations     argo     = ', argo
    write(drv%dia,*) ' Set uniform correlation radius   uniformL = ', uniformL
    write(drv%dia,*) ' Set anisotropy on corr radius    anisL    = ', anisL
    write(drv%dia,*) ' Add verbose on standard output   verbose  = ', verbose
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ''


  endif

  close(11)

  drv%sat_obs  = sat_obs
  drv%argo_obs = argo
  drv%uniformL = uniformL
  drv%anisL = anisL
  drv%Verbose = verbose

end subroutine def_nml
