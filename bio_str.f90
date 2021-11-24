MODULE bio_str

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
! Structure of biogeochemical covariance                               !
!                                                                      !
! Version 1: A.Teruzzi 2012                                            !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE bio_t

        REAL(r8),    POINTER  ::  cquot(:,:,:,:,:)     ! Component quotas
        REAL(r8),    POINTER  ::  pquot(:,:,:,:)       ! Phytoplankton component quotas
        REAL(r8),    POINTER  ::  phy(:,:,:,:,:)       ! biogeochemical variables 
        REAL(r8),    POINTER  ::  phy_ad(:,:,:,:,:)    ! biogeochemical adjoint variables
        REAL(r8),    POINTER  ::  InitialChl(:,:,:)    ! initial amount of chlorophyll
        REAL(r8),    POINTER  ::  InitialNut(:,:,:,:)  ! initial amount of nutrients
        REAL(r8),    POINTER  ::  covn3n_n1p(:,:,:)    ! covariance n3n n1p
        REAL(r8),    POINTER  ::  covn3n_chl(:,:,:)    ! covariance n3n chl
        REAL(r8),    POINTER  ::  covn1p_chl(:,:,:)    ! covariance n3n chl

        INTEGER               ::  nphy                 ! number of phytoplankton types
        INTEGER               ::  ncmp                 ! No. of phytoplankton components

        LOGICAL               ::  ApplyConditions      ! Apply conditions in snutell operations

        INTEGER(i4)           ::  N3n                  ! N3n assimilation
        INTEGER(i4)           ::  updateN1p            ! N1p update based on N3n assimilation
        INTEGER(i4)           ::  O2o                  ! O2o assimilation

   END TYPE bio_t

   TYPE (bio_t)                 :: bio

END MODULE bio_str
