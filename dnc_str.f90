MODULE dnc_str
  
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
  ! Density increment vectors                                            !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public
  
  ! ---
  ! Density increment vector for ARGO floats
  TYPE dnc_t

     INTEGER(i8)              ::  no         ! Number of all observations
     INTEGER(i8)              ::  nc         ! Number of good observations
     INTEGER(i8)              ::  k          ! Density increment index
    !  REAL(r8)                 ::  dep        ! Minimum depth for observations
    !  INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
    !  INTEGER(i8), POINTER     ::  ino(:)     ! Float number
    !  INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (0-chl, 1-N3n, 2-O2o)
     INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
     INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
     REAL(r8),    POINTER     ::  lon(:)     ! Longitute
     REAL(r8),    POINTER     ::  lat(:)     ! Latitude
     REAL(r8),    POINTER     ::  dpt(:)     ! Depth
    !  REAL(r8),    POINTER     ::  tim(:)     ! Time
     REAL(r8),    POINTER     ::  inc(:)     ! Increments
     REAL(r8),    POINTER     ::  corr(:)    ! Correlations
     REAL(r8),    POINTER     ::  err(:)     ! Nitrate std (error)
     REAL(r8),    POINTER     ::  std(:)     ! Density std
    !  REAL(r8),    POINTER     ::  res(:)     ! residual
     INTEGER(i8), POINTER     ::  ib(:)      ! i index of the nearest west point
     REAL(r8)   , POINTER     ::  pb(:)      ! distance from the nearest west point
     INTEGER(i8), POINTER     ::  jb(:)      ! j index of the nearest south point
     REAL(r8)   , POINTER     ::  qb(:)      ! distance from the nearest south point
     INTEGER(i8), POINTER     ::  kb(:)      ! k index of the nearest point below
     REAL(r8)   , POINTER     ::  rb(:)      ! distance from the nearest point below
     REAL(r8)   , POINTER     ::  pq1(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq2(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq3(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq4(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq5(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq6(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq7(:)     ! Interpolation parameter for masked grids
     REAL(r8)   , POINTER     ::  pq8(:)     ! Interpolation parameter for masked grids

     INTEGER(i4)              ::  nc_global  ! Number of global good observations

  END TYPE dnc_t

  TYPE (dnc_t)                 :: dnc

  
END MODULE dnc_str
