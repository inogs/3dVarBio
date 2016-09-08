MODULE ctl_str
  
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
  ! Cost function, control vector and optimisation arrays                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public
  
  ! ---
  ! Structure for lbfgs
  
  TYPE lbfgs_t
     
     INTEGER(i4)               ::  n          ! size of the optimisation vector
     INTEGER(i4)               ::  m          ! number of copies to be saved
     CHARACTER(LEN=60)         ::  task, csave
     LOGICAL, DIMENSION(4)     ::  lsave
     INTEGER(i4), DIMENSION(44)::  isave
! #ifndef _USE_MPI
     INTEGER(i4), POINTER      ::  nbd(:), iwa(:)
! #endif
     INTEGER(i4)               ::  iprint       
     REAL(r8)                  ::  f_b        ! The background cost function
     REAL(r8)                  ::  f_o        ! The observational cost function
     real(r8)          ::  f_c, factr ! The cost function, accuracy
     real(r8)          ::  pgtol, pgper ! Stopping criteria, percentage of initial gradient
     real(r8),  &
          DIMENSION(29)    ::  dsave
     real(r8),  &
          POINTER          ::  x_c(:)     ! The control vector (background - analyses)
     real(r8),  &
          POINTER          ::  g_c(:)     ! The gradient of f_c 
     real(r8),  &
          POINTER          ::  l_c(:), u_c(:)
! #ifndef _USE_MPI
     real(r8),  &
          POINTER          ::  wa(:), ws(:,:), wy(:,:), sy(:,:), &
          ss(:,:), wt(:,:), wn(:,:), snd(:,:),      &
          z_c(:), r_c(:), d_c(:), t_c(:)        ! Working arrays 
! #endif
     INTEGER(i4)                :: n_global   ! global size of the optimization vector (n is the local size)
     
  END TYPE lbfgs_t
  
  TYPE (lbfgs_t)                 :: ctl
  
END MODULE ctl_str
