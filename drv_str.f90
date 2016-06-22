MODULE drv_str

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
! Structure for the driver of the outer loop                           !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE drv_t

        ! CHARACTER(LEN=12)     ::  flag         ! Flag for the analysis
        INTEGER(i4)           ::  dia          ! No. of diagnostic output file
        INTEGER(i4)           ::  ntr          ! No. of outer iterations 
        INTEGER(i4)           ::  ktr          ! Outer iteration
        INTEGER(i4)           ::  im           ! Dimension of the coarse grid
        INTEGER(i4)           ::  jm           ! Dimension of the coarse grid
        INTEGER(i4), POINTER  ::  grid(:)      ! grid number for the current iterration
        LOGICAL    , POINTER  ::  read_grd(:)  ! Flag to read the grid
        REAL(r8),    POINTER  ::  ratco(:)     ! Ratio between model grid and the current grid
        REAL(r8),    POINTER  ::  ratio(:)     ! Ratio between successive grids
        INTEGER(i4), POINTER  ::  mask(:)      ! Mask used for horizontal covariances
        ! INTEGER(i4), POINTER  ::  bmd(:)       ! 1 - run barotropic model, else - do not run
        INTEGER(i4), POINTER  ::  dda(:)       ! 1 - divergence damping in analysis, else no filter
        INTEGER(i4), POINTER  ::  ddi(:)       ! 1 - divergence damping in initialisation, else no filter
        REAL(r8)              ::  f_ci         ! Inital cost function
        REAL(r8),    POINTER  ::  ro(:,:,:)    ! Vector v
        REAL(r8)              ::  f_c          ! Cost function
        REAL(r8),    POINTER  ::  ro_ad(:,:,:) ! Observational part of the cost function gradient
        REAL(r8),    POINTER  ::    msk(:,:)   ! Mask of the old grid
        ! INTEGER(i4)           ::  biol         ! Flag for the assimilation of biological parameters
        ! INTEGER(i4)           ::  bphy         ! Flag for the assimilation of biological+physical parameters

   END TYPE drv_t

   TYPE (drv_t)              :: drv

END MODULE drv_str
