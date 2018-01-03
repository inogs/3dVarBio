MODULE eof_str
  
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
  ! Structure of EOFs                                                    !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public

  TYPE eof_t
     
     LOGICAL               ::  read_eof     ! Read EOFs from file
     INTEGER(i4)           ::  neof         ! No. of EOFs
     INTEGER(i4)           ::  nreg         ! No. of regions
     INTEGER(i4)           ::  kmt          ! No. of levels of EOFs
     REAL(r8),    POINTER  ::  evcr(:,:,:)  ! Eigenvectors on regions
     REAL(r8),    POINTER  ::  evar(:,:)    ! Eigenvalues on regions
     REAL(r8),    POINTER  ::  corr(:,:,:)  ! Corelations on regions
#ifdef opt_huge_memory
     REAL(r8),    POINTER  ::  evc(:,:,:,:) ! Eigenvectors
     REAL(r8),    POINTER  ::  eva(:,:,:)   ! Eigenvalues
#else
     REAL(r8),    POINTER  ::  evc(:,:,:)   ! Eigenvectors
     REAL(r8),    POINTER  ::  eva(:,:)     ! Eigenvalues
#endif
     
     
  END TYPE eof_t
  
  TYPE (eof_t)                 :: ros
  
END MODULE eof_str
