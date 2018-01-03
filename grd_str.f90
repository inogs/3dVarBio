MODULE grd_str

  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2018 Anna Teruzzi, OGS, Trieste                         !
  !                                                                          !
  !    This file is part of 3DVarBio.
  !    3DVarBio is based on OceanVar (Dobricic, 2006)                                          !
  !                                                                          !
  !    3DVarBio is  free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    3DVarBio is  distributed in the hope that it will be useful,           !
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
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------

  use set_knd

  implicit none

  public

  TYPE grid_t

     INTEGER(i4)              ::  grd_mod      ! Grid model

     INTEGER(i4)              ::  im           ! No. points in x direction
     INTEGER(i4)              ::  jm           ! No. points in y direction
     INTEGER(i4)              ::  km           ! No. points in z direction

     REAL(r8),    POINTER     ::  ro(:,:,:)    ! Reduced order control vector
     INTEGER(i4), POINTER     ::  reg(:,:)     ! Mask for EOF regions
     INTEGER(KIND=1),    POINTER     ::  msk(:,:,:)          ! Local Sea-Land mask for scalar points
     INTEGER(KIND=1),    POINTER     ::  global_msk(:,:,:)   ! Global Sea-land mask

     REAL(r8),    POINTER     ::  ro_ad(:,:,:)    ! Reduced order control vector adjoint

     REAL(r8),    POINTER     ::  chl(:,:,:)    ! chlorophyll
     REAL(r8),    POINTER     ::  chl_ad(:,:,:) ! chlorophyll adjoint variable

     REAL(r8),    POINTER     ::  dep(:)       ! Depth

     REAL(r8),    POINTER     ::  dx(:,:)      ! dx
     REAL(r8),    POINTER     ::  dy(:,:)      ! dy
     
     REAL(r8),    POINTER     ::  alx(:,:,:)   ! Coefficient for the positive direction of the recursive filter
     REAL(r8),    POINTER     ::  aly(:,:,:)   ! Coefficient for the positive direction of the recursive filter
     REAL(r8),    POINTER     ::  btx(:,:,:)   ! Coefficient for the negative direction of the recursive filter
     REAL(r8),    POINTER     ::  bty(:,:,:)   ! Coefficient for the negative direction of the recursive filter
     REAL(r8),    POINTER     ::  scx(:,:,:)   ! Scaling factor for x direction  !laura
     REAL(r8),    POINTER     ::  scy(:,:,:)   ! Scaling factor for y direction  !laura

     INTEGER(i4)              ::  imax         ! Maximum number of extended points
     INTEGER(i4)              ::  jmax         ! Maximum number of extended points
     INTEGER(i4), POINTER     ::  imx(:)       ! Max. no. of extended pnts at each level
     INTEGER(i4), POINTER     ::  jmx(:)       ! Max. no. of extended pnts at each level
     INTEGER(i4), POINTER     ::  istp(:,:,:)    ! Extended points
     INTEGER(i4), POINTER     ::  jstp(:,:,:)    ! Extended points
     INTEGER(i4), POINTER     ::  inx(:,:,:)   ! Pointer for extended grid
     INTEGER(i4), POINTER     ::  jnx(:,:,:)   ! Pointer for extended grid
     REAL(r8),    POINTER     ::  aex(:,:,:)   ! Alpha x direction on extended grid
     REAL(r8),    POINTER     ::  aey(:,:,:)   ! Alpha y direction on extended grid
     REAL(r8),    POINTER     ::  bex(:,:,:)   ! Beta x direction on extended grid
     REAL(r8),    POINTER     ::  bey(:,:,:)   ! Beta y direction on extended grid

     REAL(r8),    POINTER     ::  lon(:,:)       ! Longitude
     REAL(r8),    POINTER     ::  lat(:,:)       ! Latitude
     REAL(r8)                 ::  NextLongitude  ! longitude value of the process below MyId

  END TYPE grid_t

  TYPE (grid_t)                 :: grd
  INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: SurfaceWaterPoints
  INTEGER(i4) :: nSurfaceWaterPoints

END MODULE grd_str
