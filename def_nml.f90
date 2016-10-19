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

#ifdef _USE_MPI
  use mpi_str
#endif

  implicit none

  INTEGER(i4), PARAMETER    :: ngrids = 3

  LOGICAL       :: read_eof
  INTEGER(i4)   :: neof, nreg, rcf_ntr, ntr
  INTEGER(i4)   :: ctl_m
  INTEGER(i4)   :: obs_chl
  INTEGER(i4)   :: obs_vdr, bmd_ncnt
  INTEGER(i4)   :: biol, bphy, nchl
  REAL(r8)      :: rcf_L, ctl_tol, ctl_per, bmd_fc1, bmd_fc2, rcf_efc, chl_dep
  INTEGER(i4)   :: grid (ngrids)
  REAL(r8)      :: ratio(ngrids)
  INTEGER(i4)   :: mask (ngrids)
  INTEGER(i4)   :: barmd(ngrids)
  INTEGER(i4)   :: divda(ngrids)
  INTEGER(i4)   :: divdi(ngrids)
  LOGICAL       :: read_grd(ngrids)


  NAMELIST /grdlst/ ntr, grid, read_grd, ratio, mask, barmd, divda, divdi
  NAMELIST /ctllst/ ctl_m, ctl_tol, ctl_per
  NAMELIST /covlst/ neof, nreg, read_eof, rcf_ntr, rcf_L, rcf_efc
  NAMELIST /biolst/ biol, bphy, nchl, chl_dep


! -------------------------------------------------------------------
! Open a formatted file for the diagnostics
! ---

  drv%dia = 12

#ifdef _USE_MPI
  if(MyRank .eq. 0) then
#endif

  open ( drv%dia, file='OceanVar.diagnostics', form='formatted' )

#ifdef _USE_MPI
endif
#endif

!---------------------------------------------------------------------
! Open the namelist
! ---

  open(11,file='var_3d_nml',form='formatted')

  ! ---
  read(11,grdlst)

#ifdef _USE_MPI
  if(MyRank .eq. 0) then
#endif

  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) '  '
  write(drv%dia,*) '                      NAMELISTS: '
  write(drv%dia,*) '  '
  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) ' GRID NAMELIST INPUT: '
  write(drv%dia,*) ' Multigrid iterrations:                  ntr    = ', ntr
  write(drv%dia,*) ' Grids:                                 grid    = ', grid (1:ntr)
  write(drv%dia,*) ' Read grids from a file:               read_grd = ', read_grd
  write(drv%dia,*) ' Ratio:                                ratio    = ', ratio(1:ntr)
  write(drv%dia,*) ' Masks:                                 mask    = ',  mask(1:ntr)
  write(drv%dia,*) ' Run barotropic model:                 barmd    = ', barmd(1:ntr)
  write(drv%dia,*) ' Divergence damping in analysis:       divda    = ', divda(1:ntr)
  write(drv%dia,*) ' Divergence damping in initialisation: divdi    = ', divdi(1:ntr)

#ifdef _USE_MPI
endif
#endif

  drv%ntr = ntr
  ALLOCATE( drv%grid (drv%ntr))      ; drv%grid (1:drv%ntr)    = grid (1:drv%ntr)
  ALLOCATE( drv%ratco(drv%ntr))      ; drv%ratco(1:drv%ntr)    = ratio(1:drv%ntr)
  ALLOCATE( drv%ratio(drv%ntr))      ; drv%ratio               = huge(drv%ratio(1))
  ALLOCATE( drv%mask (drv%ntr))      ; drv%mask (1:drv%ntr)    = mask (1:drv%ntr)
  ALLOCATE( drv%dda(drv%ntr))        ; drv%dda  (1:drv%ntr)    = divda(1:drv%ntr)
  ALLOCATE( drv%ddi(drv%ntr))        ; drv%ddi  (1:drv%ntr)    = divdi(1:drv%ntr)

  drv%ratio(        1)    = 1.0
  if(drv%ntr.gt.1) drv%ratio(2:drv%ntr)    = drv%ratco(1:drv%ntr-1) / drv%ratco(2:drv%ntr)

! ---
  read(11,ctllst)

#ifdef _USE_MPI
  if(MyRank .eq. 0) then
#endif

  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) ' MINIMIZER NAMELIST INPUT: '
  write(drv%dia,*) ' Number of saved vectors:         ctl_m    = ', ctl_m
  write(drv%dia,*) ' Minimum gradient of J:           ctl_tol  = ', ctl_tol
  write(drv%dia,*) ' Percentage of initial gradient:  ctl_per  = ', ctl_per

#ifdef _USE_MPI
endif
#endif

       ctl%m     = ctl_m
       ctl%pgtol = ctl_tol
       ctl%pgper = ctl_per

! ---
  read(11,covlst)

#ifdef _USE_MPI
  if(MyRank .eq. 0) then
#endif

  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) ' COVARIANCE NAMELIST INPUT: '
  write(drv%dia,*) ' Number of EOFs:                  neof     = ', neof
  write(drv%dia,*) ' Number of regions:               nreg     = ', nreg
  write(drv%dia,*) ' Read EOFs from a file:           read_eof = ', read_eof
  write(drv%dia,*) ' Half number of iterations:       rcf_ntr  = ', rcf_ntr
  write(drv%dia,*) ' Horizontal correlation radius:   rcf_L    = ', rcf_L
  write(drv%dia,*) ' Extension factor for coastlines: rcf_efc  = ', rcf_efc

#ifdef _USE_MPI
endif
#endif

       ros%neof     = neof
       ros%nreg     = nreg
       ros%read_eof = read_eof
       rcf%ntr      = rcf_ntr
       rcf%L        = rcf_L
       rcf%efc      = rcf_efc

! ---
  read(11,biolst)

#ifdef _USE_MPI
  if(MyRank .eq. 0) then
#endif

  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) ' BIOLOGY NAMELIST INPUT: '
  write(drv%dia,*) ' Biological assimilation          biol     = ', biol
  write(drv%dia,*) ' Biological+physical assimilation bphy     = ', bphy
  write(drv%dia,*) ' Number of phytoplankton species  nchl     = ', nchl
  write(drv%dia,*) ' Minimum depth for chlorophyll    chl_dep  = ', chl_dep

  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) '------------------------------------------------------------'
  write(drv%dia,*) ''



  close(11)
#ifdef _USE_MPI
endif
#endif

  grd%nchl = nchl
  chl%dep  = chl_dep
  drv%argo = 0 !1
  drv%ReadDomDec = 1

end subroutine def_nml
