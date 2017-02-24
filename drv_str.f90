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
     
     INTEGER(i4)           ::  dia          ! No. of diagnostic output file
     INTEGER(i4)           ::  ntr          ! No. of outer iterations 
     INTEGER(i4)           ::  ktr          ! Outer iteration
     INTEGER(i4)           ::  im           ! Dimension of the coarse grid
     INTEGER(i4)           ::  jm           ! Dimension of the coarse grid
     INTEGER(i4), POINTER  ::  grid(:)      ! grid number for the current iterration
     REAL(r8),    POINTER  ::  ratco(:)     ! Ratio between model grid and the current grid
     REAL(r8),    POINTER  ::  ratio(:)     ! Ratio between successive grids
     INTEGER(i4), POINTER  ::  mask(:)      ! Mask used for horizontal covariances
     INTEGER(i4), POINTER  ::  dda(:)       ! 1 - divergence damping in analysis, else no filter
     INTEGER(i4), POINTER  ::  ddi(:)       ! 1 - divergence damping in initialisation, else no filter
     REAL(r8),    POINTER  ::  ro(:,:,:)    ! Vector v
     REAL(r8),    POINTER  ::  msk(:,:)     ! Mask of the old grid
     
     INTEGER(i4)           ::  MyCounter    ! Number of iteration done by Tao solver
     INTEGER(i4)           ::  sat          ! Flag for the assimilation of the satellite observations
     INTEGER(i4)           ::  argo         ! Flag for the assimilation of the argo float observations
     INTEGER(i4)           ::  ReadDomDec   ! Flag for reading the Dom_Dec_jp*.ascii files
     INTEGER(i4)           ::  Verbose      ! Flag for printing verbose output

  END TYPE drv_t
  
  TYPE (drv_t)              :: drv
  
END MODULE drv_str
