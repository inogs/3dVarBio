MODULE ctl_str
  
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
  ! Cost function, control vector and optimisation arrays                !
  !                                                                      !
  ! Version 1: A. Teruzzi 2018                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  
  implicit none
  
  public
  
  ! ---
  ! Structure for lbfgs
  
  TYPE lbfgs_t
     
     INTEGER(i4)               ::  n          ! size of the optimisation vector

     REAL(r8)                  ::  f_b        ! The background cost function
     REAL(r8)                  ::  f_o        ! The observational cost function
     real(r8)           ::  f_c, factr ! The cost function, accuracy
     real(r8)           ::  pgtol, pgper ! Stopping criteria, percentage of initial gradient

     real(r8), POINTER  ::  x_c(:)     ! The control vector (background - analyses)
     real(r8), POINTER  ::  g_c(:)     ! The gradient of f_c 
     INTEGER(i4)        :: n_global   ! global size of the optimization vector (n is the local size)
     
  END TYPE lbfgs_t
  
  TYPE (lbfgs_t)                 :: ctl
  
END MODULE ctl_str
