MODULE tao_str

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
#include "petsc/finclude/petscvecdef.h"
 use set_knd
 use petscvec
 implicit none
 
 public


! #include "tao_minimizer.h"
 ! ---
 ! Structure for lbfgs

 TYPE petsc_str
    
    INTEGER(i4)               ::  n          ! local size of the optimisation vector
    INTEGER(i4)               ::  n_global   ! global size of the optimization vector
    real(r8)                 ::  f_b        ! The background cost function
    real(r8)                 ::  f_o        ! The observational cost function
    real(r8)                 ::  f_c, factr ! The cost function, accuracy
    ! PetscReal                 ::  f_b        ! The background cost function
    ! PetscReal                 ::  f_o        ! The observational cost function
    ! PetscReal                 ::  f_c, factr ! The cost function, accuracy

    ! PetscFortranAddr                       ::  x_c        ! The control vector (background - analyses)
    ! PetscFortranAddr                       ::  g_c        ! The gradient of f_c 

    Vec                       ::  x_c        ! The control vector (background - analyses)
    Vec                       ::  g_c        ! The gradient of f_c 
    
 END TYPE petsc_str
 
 ! TYPE (petsc_str)                 :: NewCtl
 
END MODULE tao_str
