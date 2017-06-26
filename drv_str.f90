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
  ! Structure for the driver and assimilation flags                      !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public
  
  TYPE drv_t
     
     INTEGER(i4)           ::  dia          ! No. of diagnostic output file
     
     INTEGER(i4)           ::  MyCounter    ! Number of iteration done by Tao solver
     INTEGER(i4)           ::  sat_obs      ! Flag for the assimilation of the satellite observations
     INTEGER(i4)           ::  argo_obs     ! Flag for the assimilation of the argo float observations
     INTEGER(i4)           ::  chl          ! Flag for the chlorophyll assimilation
     INTEGER(i4)           ::  uniformL     ! Flag for setting uniform correlation radius (1 = non uniform)
     INTEGER(i4)           ::  anisL        ! Flag for setting anisotropy on correlation radius (1 = anisotropy)
     INTEGER(i4)           ::  Verbose      ! Flag for printing verbose output
     INTEGER(i4)           ::  nut          ! Flag for the chlorophyll assimilation

  END TYPE drv_t
  
  TYPE (drv_t)              :: drv
  
END MODULE drv_str
