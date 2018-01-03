MODULE bio_str

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
! Structure of biogeochemical covariance                               !
!                                                                      !
! Version 1: A. Teruzzi 2018                                            !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE bio_t

        REAL(r8),    POINTER  ::  cquot(:,:,:,:,:)     ! Component quotas
        REAL(r8),    POINTER  ::  pquot(:,:,:,:)       ! Phytoplankton component quotas
        REAL(r8),    POINTER  ::  phy(:,:,:,:,:)       ! biogeochemical variables 
        REAL(r8),    POINTER  ::  phy_ad(:,:,:,:,:)    ! biogeochemical adjoint variables
        REAL(r8),    POINTER  ::  InitialChl(:,:,:)     ! initial amount of chlorophyll

        INTEGER               ::  nphy                 ! number of phytoplankton types
        INTEGER               ::  ncmp                 ! No. of phytoplankton components

        LOGICAL               ::  ApplyConditions       ! Apply conditions in snutell operations
        
   END TYPE bio_t

   TYPE (bio_t)                 :: bio

END MODULE bio_str
