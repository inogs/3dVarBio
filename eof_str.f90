MODULE eof_str
  
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
  ! Structure of EOFs                                                    !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public

  TYPE eof_t
     
     LOGICAL               ::  read_eof     ! Read EOFs from file
     INTEGER(i4)           ::  neof         ! Total No. of EOFs
     INTEGER(i4)           ::  neof_chl     ! No. of EOFs for chl
     INTEGER(i4)           ::  neof_nut     ! No. of EOFs for nutrients
     INTEGER(i4)           ::  neof_n3n     ! No. of EOFs for N3n
     INTEGER(i4)           ::  neof_o2o     ! No. of EOFs for O2o
     INTEGER(i4)           ::  neof_multi   ! No. of EOFs for multivariate
     INTEGER(i4)           ::  nreg         ! No. of regions
     INTEGER(i4)           ::  kmt          ! No. of levels of EOFs
     INTEGER(i4)           ::  kmchl  ! No. of levels of multi EOFs for chl
     INTEGER(i4)           ::  kmnit  ! No. of levels of multi EOFs for nit
     REAL(r8),    POINTER  ::  evcr(:,:,:)  ! Eigenvectors on regions
     REAL(r8),    POINTER  ::  evar(:,:)    ! Eigenvalues on regions
     REAL(r8),    POINTER  ::  corr(:,:,:)  ! Corelations on regions
#ifdef opt_huge_memory
     REAL(r8),    POINTER  ::  evc(:,:,:,:) ! Eigenvectors
     REAL(r8),    POINTER  ::  eva(:,:,:)   ! Eigenvalues
#else
     REAL(r8),    POINTER  ::  evc_chl(:,:,:)   ! Eigenvectors
     REAL(r8),    POINTER  ::  eva_chl(:,:)     ! Eigenvalues
     REAL(r8),    POINTER  ::  evc_n3n(:,:,:)   ! Eigenvectors
     REAL(r8),    POINTER  ::  eva_n3n(:,:)     ! Eigenvalues
     REAL(r8),    POINTER  ::  evc_o2o(:,:,:)   ! Eigenvectors
     REAL(r8),    POINTER  ::  eva_o2o(:,:)     ! Eigenvalues
     REAL(r8),    POINTER  ::  evc_multi(:,:,:)   ! Eigenvectors
     REAL(r8),    POINTER  ::  eva_multi(:,:)     ! Eigenvalues
#endif
     
     
  END TYPE eof_t
  
  TYPE (eof_t)                 :: ros
  
END MODULE eof_str
