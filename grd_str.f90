MODULE grd_str

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
! Structure of the grid                                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE grid_t

        INTEGER(i4)              ::  grd_mod      ! Grid model
        LOGICAL                  ::  read_grd     ! Read the grid from a file

        INTEGER(i4)              ::  im           ! No. points in x direction
        INTEGER(i4)              ::  jm           ! No. points in y direction
        INTEGER(i4)              ::  km           ! No. points in z direction
        INTEGER(i4)              ::  nps          ! No. of ocean points 

        REAL(r8)                 ::  dln          ! Resolution in the x direction
        REAL(r8)                 ::  dlt          ! Resolution in the y direction

        REAL(r8),    POINTER     ::  ro(:,:,:)    ! Reduced order control vector
        INTEGER(i4), POINTER     ::  reg(:,:)     ! Mask for EOF regions
        REAL(r8),    POINTER     ::  msk(:,:,:)   ! Sea-Land mask for scalar points
        REAL(r8),    POINTER     ::  ums(:,:,:)   ! Sea-Land mask for u points
        REAL(r8),    POINTER     ::  vms(:,:,:)   ! Sea-Land mask for v points
        REAL(r8),    POINTER     ::  hgt(:,:)     ! Topography
        REAL(r8),    POINTER     ::    f(:,:)     ! Coriolis term

        ! REAL(r8),    POINTER     ::  tem(:,:,:)   ! Temperature increment
        ! REAL(r8),    POINTER     ::  sal(:,:,:)   ! Salinity increment
        ! REAL(r8),    POINTER     ::  uvl(:,:,:)   ! u componnet of velocity increment
        ! REAL(r8),    POINTER     ::  vvl(:,:,:)   ! v componnet of velocity increment
        ! REAL(r8),    POINTER     ::  eta(:,:)     ! Sea level increment

        ! REAL(r8),    POINTER     ::  temb(:,:,:)   ! Temperature background
        ! REAL(r8),    POINTER     ::  salb(:,:,:)   ! Salinity background
        REAL(r8),    POINTER     ::  uvlb(:,:,:)   ! u componnet of velocity background
        REAL(r8),    POINTER     ::  vvlb(:,:,:)   ! v componnet of velocity background
        ! REAL(r8),    POINTER     ::  etab(:,:)     ! Sea level background

        REAL(r8),    POINTER     ::  mdt(:,:)      ! Mean dynamic topography
        ! REAL(r8),    POINTER     ::  sla(:,:)      ! Sea level anomaly

        REAL(r8),    POINTER     ::  ro_ad(:,:,:)    ! Reduced order control vector adjoint
        ! REAL(r8),    POINTER     ::  tem_ad(:,:,:)   ! Temperature adjoint
        ! REAL(r8),    POINTER     ::  sal_ad(:,:,:)   ! Salinity adjoint
        ! REAL(r8),    POINTER     ::  uvl_ad(:,:,:)   ! u componnet of velocity adjoint
        ! REAL(r8),    POINTER     ::  vvl_ad(:,:,:)   ! v componnet of velocity adjoint
        ! REAL(r8),    POINTER     ::  eta_ad(:,:)     ! Sea level adjoint

        ! REAL(r8),    POINTER     ::  dns(:,:,:)      ! density
        ! REAL(r8),    POINTER     ::  b_x(:,:,:)      ! bouyancy force 
        ! REAL(r8),    POINTER     ::  b_y(:,:,:)      ! bouyancy force 
        REAL(r8),    POINTER     ::  bx(:,:)         ! bouyancy force integral
        REAL(r8),    POINTER     ::  by(:,:)         ! bouyancy force integral

        INTEGER(i4)              ::  nchl            ! No. of phytoplankton species
        REAL(r8),    POINTER     ::  chl(:,:,:,:)    ! chlorophyll
        REAL(r8),    POINTER     ::  chl_ad(:,:,:,:) ! chlorophyll adjoint variable

        REAL(r8),    POINTER     ::  lon(:,:)       ! Longitude
        REAL(r8),    POINTER     ::  lat(:,:)       ! Latitude
        REAL(r8),    POINTER     ::  dep(:)       ! Depth

        REAL(r8),    POINTER     ::  hvst(:,:,:)  ! Horizontal diff. coef. for temperature 
        REAL(r8),    POINTER     ::  hvss(:,:,:)  ! Horizontal diff. coef. for salinity   
        REAL(r8),    POINTER     ::  hvsp(:,:,:)  ! Horizontal diff. coef. for stream function
        REAL(r8),    POINTER     ::  vvst(:,:,:)  ! Vertical diff. coef. for temperature 
        REAL(r8),    POINTER     ::  vvss(:,:,:)  ! Vertical diff. coef. for salinity   
        REAL(r8),    POINTER     ::  vvsp(:,:,:)  ! Vertical diff. coef. for stream function
        REAL(r8),    POINTER     ::  dx(:,:)      ! dx
        REAL(r8),    POINTER     ::  dy(:,:)      ! dy
        REAL(r8),    POINTER     ::  dz(:)        ! dz
        REAL(r8),    POINTER     ::  dxdy(:,:)    ! dx*dy
        REAL(r8)                 ::  adxdy        ! Mean dx*dy


        REAL(r8),    POINTER     ::  alx(:,:)     ! Coefficient for the positive direction of the recursive filter
        REAL(r8),    POINTER     ::  aly(:,:)     ! Coefficient for the positive direction of the recursive filter
        REAL(r8),    POINTER     ::  btx(:,:)     ! Coefficient for the negative direction of the recursive filter
        REAL(r8),    POINTER     ::  bty(:,:)     ! Coefficient for the negative direction of the recursive filter
        REAL(r8),    POINTER     ::  scx(:,:)     ! Scaling factor for x direction
        REAL(r8),    POINTER     ::  scy(:,:)     ! Scaling factor for y direction
        REAL(r8),    POINTER     ::  msr(:,:,:)   ! Sea-land mask used in the recursive filter
        INTEGER(i4)              ::  imax         ! Maximum number of extended points
        INTEGER(i4)              ::  jmax         ! Maximum number of extended points
        INTEGER(i4), POINTER     ::  imx(:)       ! Max. no. of extended pnts at each level
        INTEGER(i4), POINTER     ::  jmx(:)       ! Max. no. of extended pnts at each level
        INTEGER(i4), POINTER     ::  istp(:,:)    ! Extended points
        INTEGER(i4), POINTER     ::  jstp(:,:)    ! Extended points
        INTEGER(i4), POINTER     ::  inx(:,:,:)   ! Pointer for extended grid
        INTEGER(i4), POINTER     ::  jnx(:,:,:)   ! Pointer for extended grid
        REAL(r8),    POINTER     ::  fct(:,:,:)   ! Normalisation factor
        REAL(r8),    POINTER     ::  aex(:,:,:)   ! Alpha x direction on extended grid
        REAL(r8),    POINTER     ::  aey(:,:,:)   ! Alpha y direction on extended grid
        REAL(r8),    POINTER     ::  bex(:,:,:)   ! Beta x direction on extended grid
        REAL(r8),    POINTER     ::  bey(:,:,:)   ! Beta y direction on extended grid


   END TYPE grid_t

   TYPE (grid_t)                 :: grd
   INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: SurfaceWaterpoints
   INTEGER(i4) :: nSurfaceWaterPoints
   REAL(r4),   ALLOCATABLE, DIMENSION(:,:,:)     ::  Dump_chl, Dump_vip
   REAL(r4),   ALLOCATABLE, DIMENSION(:,:)       ::  Dump_msk
END MODULE grd_str
