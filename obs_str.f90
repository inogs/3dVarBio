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
        INTEGER(i8)              ::  sla        ! Flag for assimilation of SLA
        INTEGER(i8)              ::  arg        ! Flag for assimilation of ARGO floats
        INTEGER(i8)              ::  xbt        ! Flag for assimilation of XBTs
        INTEGER(i8)              ::  gld        ! Flag for assimilation of gliders
        INTEGER(i8)              ::  tra        ! Flag for assimilation of Argo trajectories
        INTEGER(i8)              ::  trd        ! Flag for assimilation of driftres trajectories
        INTEGER(i8)              ::  vdr        ! Flag for assimilation of velocities by drifters
        INTEGER(i8)              ::  gvl        ! Flag for assimilation of velocities by gliders
        INTEGER(i8)              ::  chl        ! Flag for assimilation of chlorophyll

   END TYPE obs_t

   TYPE (obs_t)                 :: obs

! ---
! Observational vector for SLA
   TYPE sla_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
        INTEGER(i8), POINTER     ::  ino(:)     ! Instrument
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

   END TYPE sla_t

   TYPE (sla_t)                 :: sla

! ---
! Observational vector for ARGO floats
   TYPE arg_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (1-temperature, 2-salinity)
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        REAL(r8),    POINTER     ::  lon(:)     ! Longitute
        REAL(r8),    POINTER     ::  lat(:)     ! Latitude
        REAL(r8),    POINTER     ::  dpt(:)     ! Depth
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

   END TYPE arg_t

   TYPE (arg_t)                 :: arg

! ---
! Observational vector for XBT profiles
   TYPE xbt_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (1-temperature)
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        REAL(r8),    POINTER     ::  lon(:)     ! Longitute
        REAL(r8),    POINTER     ::  lat(:)     ! Latitude
        REAL(r8),    POINTER     ::  dpt(:)     ! Depth
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

   END TYPE xbt_t

   TYPE (xbt_t)                 :: xbt

! ---
! Observational vector for gliders
   TYPE gld_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8)              ::  kdp        ! Model level corresponding to dep
        INTEGER(i8), POINTER     ::  ino(:)     ! Glider number
        INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (1-temperature, 2-salinity)
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        REAL(r8),    POINTER     ::  lon(:)     ! Longitute
        REAL(r8),    POINTER     ::  lat(:)     ! Latitude
        REAL(r8),    POINTER     ::  dpt(:)     ! Depth
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

   END TYPE gld_t

   TYPE (gld_t)                 :: gld

! ---
! Observational vector for Argo trajectories
   TYPE tra_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        INTEGER(i8)              ::  nt         ! Number of time steps
        INTEGER(I8)              ::  im         ! I dimension of the grid
        INTEGER(I8)              ::  jm         ! J dimension of the grid
        INTEGER(I8)              ::  km         ! K dimension of the grid
        REAL(r8)                 ::  dpt        ! Depth of observations
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        REAL(r8),    POINTER     ::  dtm(:)     ! Time spent under surface
        REAL(r8),    POINTER     ::  loi(:)     ! Initial observed longitude
        REAL(r8),    POINTER     ::  lai(:)     ! Initial observed latitude
        REAL(r8),    POINTER     ::  lof(:)     ! Final observed longitude
        REAL(r8),    POINTER     ::  laf(:)     ! Final observed latitude
        REAL(r8),    POINTER     ::  lob(:,:)   ! Longitudes of simulated positions
        REAL(r8),    POINTER     ::  lab(:,:)   ! Latitudes of simulated positions
        REAL(r8),    POINTER     ::  loa(:)     ! Longitudes of analysed positions
        REAL(r8),    POINTER     ::  laa(:)     ! Latitudes of analysed positions
        REAL(r8),    POINTER     ::  xob(:)     ! Observed value in x direction
        REAL(r8),    POINTER     ::  erx(:)     ! Observational error in x direction
        REAL(r8),    POINTER     ::  rex(:)     ! Residual in x direction
        REAL(r8),    POINTER     ::  inx(:)     ! Increments in x direction
        REAL(r8),    POINTER     ::  yob(:)     ! Observed value in y direction
        REAL(r8),    POINTER     ::  ery(:)     ! Observational error in y direction
        REAL(r8),    POINTER     ::  rey(:)     ! Residual in y direction
        REAL(r8),    POINTER     ::  iny(:)     ! Increments in y direction
        REAL(r8),    POINTER     ::  err(:)     ! Observational error in meters
        REAL(r8)   , POINTER     ::  umn(:,:)   ! U component of background velocity
        REAL(r8)   , POINTER     ::  vmn(:,:)   ! V component of background velocity
        REAL(r8)   , POINTER     ::  dx(:,:)    ! Delta x on trajectory grid
        REAL(r8)   , POINTER     ::  dy(:,:)    ! Delta x on trajectory grid
        REAL(r8)   , POINTER     ::  lon(:,:)   ! Longitudes of grid points
        REAL(r8)   , POINTER     ::  lat(:,:)   ! Latitudes of grid points
        REAL(r8)   , POINTER     ::  xmn(:,:)   ! X coordinate of mean trajectory position
        REAL(r8)   , POINTER     ::  ymn(:,:)   ! Y coordinate of mean trajectory position
        REAL(r8)   , POINTER     ::  tim(:)     ! Time of the duration of the trajectory
        REAL(r8)   , POINTER     ::  xtl(:)     ! Delta X of trajectory correction
        REAL(r8)   , POINTER     ::  ytl(:)     ! Delta Y of trajectory correction
        REAL(r8)   , POINTER     ::  xtl_ad(:)  ! Delta X of trajectory correction (adjoint)
        REAL(r8)   , POINTER     ::  ytl_ad(:)  ! Delta Y of trajectory correction (adjoint)
        INTEGER(i8), POINTER     ::  i1(:,:)    ! i index for interpolation between grids
        INTEGER(i8), POINTER     ::  j1(:,:)    ! j index for interpolation between grids
        REAL(r8)   , POINTER     ::  pq1(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq2(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq3(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq4(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  uvl(:,:)   ! Delta u on trajectory grid
        REAL(r8)   , POINTER     ::  vvl(:,:)   ! Delta v on trajectory grid
        REAL(r8)   , POINTER     ::  uvl_ad(:,:)! Delta u on trajectory grid (adjoint)
        REAL(r8)   , POINTER     ::  vvl_ad(:,:)! Delta v on trajectory grid (adjoint)


   END TYPE tra_t

   TYPE (tra_t)                 :: tra

! ---
! Observational vector for trajectories of surface drifters
   TYPE trd_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        INTEGER(i8)              ::  nt         ! Number of time steps
        INTEGER(I8)              ::  im         ! I dimension of the grid
        INTEGER(I8)              ::  jm         ! J dimension of the grid
        INTEGER(I8)              ::  km         ! K dimension of the grid
        REAL(r8)                 ::  dpt        ! Depth of observations
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        REAL(r8),    POINTER     ::  dtm(:)     ! Time spent under surface
        REAL(r8),    POINTER     ::  loi(:)     ! Initial observed longitude
        REAL(r8),    POINTER     ::  lai(:)     ! Initial observed latitude
        REAL(r8),    POINTER     ::  lof(:)     ! Final observed longitude
        REAL(r8),    POINTER     ::  laf(:)     ! Final observed latitude
        REAL(r8),    POINTER     ::  lob(:,:)   ! Longitudes of simulated positions
        REAL(r8),    POINTER     ::  lab(:,:)   ! Latitudes of simulated positions
        REAL(r8),    POINTER     ::  loa(:)     ! Longitudes of analysed positions
        REAL(r8),    POINTER     ::  laa(:)     ! Latitudes of analysed positions
        REAL(r8),    POINTER     ::  xob(:)     ! Observed value in x direction
        REAL(r8),    POINTER     ::  erx(:)     ! Observational error in x direction
        REAL(r8),    POINTER     ::  rex(:)     ! Residual in x direction
        REAL(r8),    POINTER     ::  inx(:)     ! Increments in x direction
        REAL(r8),    POINTER     ::  yob(:)     ! Observed value in y direction
        REAL(r8),    POINTER     ::  ery(:)     ! Observational error in y direction
        REAL(r8),    POINTER     ::  rey(:)     ! Residual in y direction
        REAL(r8),    POINTER     ::  iny(:)     ! Increments in y direction
        REAL(r8),    POINTER     ::  err(:)     ! Observational error in meters
        REAL(r8)   , POINTER     ::  umn(:,:)   ! U component of background velocity
        REAL(r8)   , POINTER     ::  vmn(:,:)   ! V component of background velocity
        REAL(r8)   , POINTER     ::  dx(:,:)    ! Delta x on trajectory grid
        REAL(r8)   , POINTER     ::  dy(:,:)    ! Delta x on trajectory grid
        REAL(r8)   , POINTER     ::  lon(:,:)   ! Longitudes of grid points
        REAL(r8)   , POINTER     ::  lat(:,:)   ! Latitudes of grid points
        REAL(r8)   , POINTER     ::  xmn(:,:)   ! X coordinate of mean trajectory position
        REAL(r8)   , POINTER     ::  ymn(:,:)   ! Y coordinate of mean trajectory position
        REAL(r8)   , POINTER     ::  tim(:)     ! Time of the duration of the trajectory
        REAL(r8)   , POINTER     ::  xtl(:)     ! Delta X of trajectory correction
        REAL(r8)   , POINTER     ::  ytl(:)     ! Delta Y of trajectory correction
        REAL(r8)   , POINTER     ::  xtl_ad(:)  ! Delta X of trajectory correction (adjoint)
        REAL(r8)   , POINTER     ::  ytl_ad(:)  ! Delta Y of trajectory correction (adjoint)
        INTEGER(i8), POINTER     ::  i1(:,:)    ! i index for interpolation between grids
        INTEGER(i8), POINTER     ::  j1(:,:)    ! j index for interpolation between grids
        REAL(r8)   , POINTER     ::  pq1(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq2(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq3(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  pq4(:,:)   ! Parameter for interpolation between grids
        REAL(r8)   , POINTER     ::  uvl(:,:)   ! Delta u on trajectory grid
        REAL(r8)   , POINTER     ::  vvl(:,:)   ! Delta v on trajectory grid
        REAL(r8)   , POINTER     ::  uvl_ad(:,:)! Delta u on trajectory grid (adjoint)
        REAL(r8)   , POINTER     ::  vvl_ad(:,:)! Delta v on trajectory grid (adjoint)


   END TYPE trd_t

   TYPE (trd_t)                 :: trd

! ---
! Observational vector for velocity from drifters
   TYPE vdr_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (1 - u component, 2 - v component)
        REAL(r8),    POINTER     ::  lon(:)     ! Longitude
        REAL(r8),    POINTER     ::  lat(:)     ! Latitude
        REAL(r8),    POINTER     ::  dpt(:)     ! Depth
        INTEGER(i8), POINTER     ::  kdp(:)     ! Model level corresponding to dep
        REAL(r8),    POINTER     ::  tim(:)     ! Time
        REAL(r8),    POINTER     ::  tms(:)     ! Starting time for averaging
        REAL(r8),    POINTER     ::  tme(:)     ! Final time for averaging
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
        INTEGER(i8), POINTER     ::  nav(:)      ! Number of time steps for averaging
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

   END TYPE vdr_t

   TYPE (vdr_t)                 :: vdr

! ---
! Observational vector for velocity from gliders
   TYPE gvl_t

        INTEGER(i8)              ::  no         ! Number of all observations
        INTEGER(i8)              ::  nc         ! Number of good observations
        REAL(r8)                 ::  dep        ! Minimum depth for observations
        INTEGER(i8), POINTER     ::  flg(:)     ! Quality flag
        INTEGER(i8), POINTER     ::  flc(:)     ! Temporary flag for multigrid
        INTEGER(i8), POINTER     ::  ino(:)     ! Float number
        INTEGER(i8), POINTER     ::  par(:)     ! Parameter flag (1 - u component, 2 - v component)
        REAL(r8),    POINTER     ::  lon(:)     ! Longitude
        REAL(r8),    POINTER     ::  lat(:)     ! Latitude
        REAL(r8),    POINTER     ::  dpt(:)     ! Depth
        REAL(r8),    POINTER     ::  dzr(:,:)   ! Relative thickness
        INTEGER(i8), POINTER     ::  kdp(:)     ! Model level corresponding to dep
        REAL(r8),    POINTER     ::  tim(:)     ! Time
        REAL(r8),    POINTER     ::  tms(:)     ! Starting time for averaging
        REAL(r8),    POINTER     ::  tme(:)     ! Final time for averaging
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
        INTEGER(i8), POINTER     ::  nav(:)      ! Number of time steps for averaging
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

   END TYPE gvl_t

   TYPE (gvl_t)                 :: gvl

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
