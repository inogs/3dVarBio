MODULE set_knd

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
! The precision of reals and integers                                  !
!                                                                      !
! Version 1: A. Teruzzi 2018                                           !
!-----------------------------------------------------------------------


implicit none

public

   INTEGER, PARAMETER ::                &
      r4 = SELECTED_REAL_KIND( 6, 37),  &  ! real*4
      r8 = SELECTED_REAL_KIND(12,307)      ! real*8

   INTEGER, PARAMETER ::                &
      i4 = SELECTED_INT_KIND(9) ,       &  ! integer*4
      i8 = i4! ELECTED_INT_KIND(14)           ! integer*8

  type DoubleGrid
    sequence
    real(r8) chl
    real(r8) chl_ad
  end type DoubleGrid

END MODULE set_knd
