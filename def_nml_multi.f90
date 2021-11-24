subroutine def_nml_multi

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2018 Anna Teruzzi, OGS, Trieste                             !
!                                                                          !
!    This file is part of 3DVarBio.                                        !
!    3DVarBio is based on OceanVar (Dobricic, 2006)                        !
!                                                                          !
!    3DVarBio is  free software: you can redistribute it and/or modify.    !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    3DVarBio is  distributed in the hope that it will be useful,          !
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
! Define analysis parameters from namelists                            !
! for multi platform and multivariate DA                               !
! (other general DA parameters are defined in def_nml.f90)             !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
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
  use da_params

  implicit none

  LOGICAL       :: ApplyConditions
  INTEGER(i4)   :: chl_assim, chl_upnut, nut, multiv, N3n, O2o, updateN1p
  INTEGER(i4)   :: nphyto, uniformL, anisL
  REAL(r8)      :: chl_dep
  INTEGER(i4)   :: argo, sat_obs, ncmp
  
  !NAMELIST /ctllst/ ctl_tol, ctl_per
  !NAMELIST /covlst/ neof_chl, neof_n3n, neof_o2o, nreg, read_eof, rcf_ntr, rcf_L, rcf_efc
  NAMELIST /biolst/ chl_assim, chl_upnut, nut, multiv, nphyto, chl_dep, ncmp, ApplyConditions, N3n, updateN1p, O2o
  NAMELIST /params/ sat_obs, argo, uniformL, anisL


! -------------------------------------------------------------------
! Open a formatted file for the diagnostics
! ---

  drv%dia = 12

  if(MyId .eq. 0) then
    open ( drv%dia, file='DA__FREQ_1/OceanVar.dia_multinml.'//DA_DATE, form='formatted' )
  endif

!---------------------------------------------------------------------
! Open the namelist
! ---

  open(11,file='DA__FREQ_1/satfloat.'//DA_DATE//'.nml',form='formatted')

  ! ---

  read(11,biolst)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' BIOLOGY NAMELIST INPUT: '
    write(drv%dia,*) ' Chlorophyll assimilation             chl_assim = ', chl_assim
    write(drv%dia,*) ' N3n update based on chl assimilation chl_upnut = ', chl_upnut
    write(drv%dia,*) ' Nutrient assimilation                      nut = ', nut
    write(drv%dia,*) ' Multivariate assimilation               multiv = ', multiv
    write(drv%dia,*) ' Number of phytoplankton species          nphyt = ', nphyto
    write(drv%dia,*) ' Minimum depth for chlorophyll          chl_dep = ', chl_dep
    write(drv%dia,*) ' Number of phytoplankton components        ncmp = ', ncmp
    write(drv%dia,*) ' Apply conditions flag          ApplyConditions = ', ApplyConditions
    write(drv%dia,*) ' N3n assimilation                           N3n = ', N3n
    write(drv%dia,*) ' N1p update based on N3n assimilation updateN1p = ', updateN1p
    write(drv%dia,*) ' O2o assimilation                           O2o = ', O2o

  endif

  drv%chl_assim  = chl_assim
  drv%chl_upnut = chl_upnut
  drv%nut  = nut
  drv%multiv = multiv
  bio%nphy = nphyto
  sat%dep  = chl_dep
  bio%ncmp = ncmp
  bio%ApplyConditions = ApplyConditions
  bio%N3n = N3n
  bio%updateN1p = updateN1p
  bio%O2o = O2o

  ! ---
  
  read(11,params)

  if(MyId .eq. 0) then

    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ' PARAMETERS NAMELIST INPUT: '
    write(drv%dia,*) ' Read Satellite observations      sat_obs  = ', sat_obs
    write(drv%dia,*) ' Read ARGO float observations     argo     = ', argo
    write(drv%dia,*) ' Set uniform correlation radius   uniformL = ', uniformL
    write(drv%dia,*) ' Set anisotropy on corr radius    anisL    = ', anisL
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) '------------------------------------------------------------'
    write(drv%dia,*) ''


  endif

  close(11)

  drv%sat_obs  = sat_obs
  drv%argo_obs = argo
  drv%uniformL = uniformL
  drv%anisL = anisL

end subroutine def_nml_multi
