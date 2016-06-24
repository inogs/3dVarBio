MODULE obs_str
  
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
  ! Observational vectors                                                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public
  
  ! ---
  ! Observational vector in the cost function
  TYPE obs_t
     
     INTEGER(i8)              ::  no         ! Number of observations
     INTEGER(i8)              ::  k          ! Observation index 
     REAL(r8),    POINTER     ::  inc(:)     ! Increments
     REAL(r8),    POINTER     ::  amo(:)     ! Analysis - observation
     REAL(r8),    POINTER     ::  res(:)     ! residual
     REAL(r8),    POINTER     ::  err(:)     ! Observational error
     REAL(r8),    POINTER     ::  gra(:)     ! Observational gradient
     
  END TYPE obs_t
  
  TYPE (obs_t)                 :: obs
  
  ! ---
  ! Observational vector for Chlorophyll
  TYPE chl_t
     
     INTEGER(i8)              ::  no         ! Number of all observations
     INTEGER(i8)              ::  nc         ! Number of good observations
     REAL(r8)                 ::  dep        ! Minimum depth for observations
     INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
     INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
     INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
     REAL(r8),    POINTER     ::  lon(:)     ! Longitute
     REAL(r8),    POINTER     ::  lat(:)     ! Latitude
     REAL(r8),    POINTER     ::  tim(:)     ! Time
     REAL(r8),    POINTER     ::  val(:)     ! Observed value
     REAL(r8),    POINTER     ::  bac(:)     ! Background value
     REAL(r8),    POINTER     ::  inc(:)     ! Increments
     REAL(r8),    POINTER     ::  bia(:)     ! Bias
     REAL(r8),    POINTER     ::  err(:)     ! Observational error
     REAL(r8),    POINTER     ::  res(:)     ! residual
     REAL(r8),    POINTER     ::  b_a(:)     ! Background - analyses
     INTEGER(i8), POINTER     ::  ib(:)      ! i index of the nearest west point
     REAL(r8)   , POINTER     ::  pb(:)      ! distance from the nearest west point
     INTEGER(i8), POINTER     ::  jb(:)      ! j index of the nearest south point
     REAL(r8)   , POINTER     ::  qb(:)      ! distance from the nearest south point
     REAL(r8)   , POINTER     ::  pq1(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq2(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq3(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq4(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  dpt(:)     ! Maximum depth of surrounding points
     INTEGER(i8), POINTER     ::  kb(:)      ! k index of bottom point for vertical integration
     REAL(r8),    POINTER     ::  dzr(:,:)   ! Relative thickness
     
  END TYPE chl_t
  
  TYPE (chl_t)                 :: chl
  
  ! ---

  
END MODULE obs_str
