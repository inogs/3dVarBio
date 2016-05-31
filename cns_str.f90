MODULE cns_str

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
! Structure of constants                                               !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE rcf_t

        INTEGER(i4)          ::  ntr     ! No. of iterations (half of)
        REAL(r8)             ::  dx      ! Grid resolution (m)
        REAL(r8)             ::  L       ! Correlation radius
        REAL(r8)             ::  E       ! Norm
        REAL(r8)             ::  alp     ! Filter weight
        INTEGER(i4)          ::  ntb     ! Number of points in the table
        REAL(r8)             ::  dsmn    ! Minimum distance 
        REAL(r8)             ::  dsmx    ! Maximum distance 
        REAL(r8)             ::  dsl     ! Table increment
        REAL(r8), POINTER    ::  al(:)   ! Filter weights in the table
        REAL(r8), POINTER    ::  sc(:)   ! Filter scaling factors in the table
        REAL(r8)             ::  scl     ! Scaling factor
        REAL(r8)             ::  efc     ! Scaling factor for extended points
        INTEGER(i4)          ::  kstp    ! Step for extended points

   END TYPE rcf_t

   TYPE (rcf_t)              :: rcf

END MODULE cns_str
